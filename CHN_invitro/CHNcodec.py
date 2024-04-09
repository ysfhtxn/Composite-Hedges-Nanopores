from base64 import decode
import numpy as np
from copy import deepcopy
from typing import List
from utils import composite_letter, Hypothesis, code_pattern, HypothesisTree
from scipy.stats import wasserstein_distance
from scipy.spatial import distance


DEBUG0 = False
DEBUG1 = False

"""
This version added GC content and max runs of homopolymer constrant.
Date : 
"""
def ran_hash(u):
    """
    哈希编码
    :param u: 编码值
    :return: 哈希值
    """
    v = u * 3935559000370003845 + 2691343689449507681
    v ^= v >> 21
    v ^= v << 37
    v ^= v >> 4
    v *= 4768777513237032717
    v ^= v << 20
    v ^= v >> 41
    v ^= v << 5
    return v

def unpack_vbits(msg, pattern):
    """
    将msg转化成二进制数据, little-endian
    :param msg: 输入数据
    :param pattern: 进行哈希变换时的数据位
    :return: 二进制数据
    """
    vbits = []
    for word in msg:
        mask = 1
        for _ in range(pattern):
            vbits.append((word & mask))
            word = word >> 1
    return vbits

def pack_vbits(vbits, pattern):
    """
    将二进制数据转化成msg, little-endian
    :param vbits: 输入数据
    :param pattern: 进行哈希变换时的数据位
    :return: 二进制数据
    """
    msg = []
    for i in range(0, len(vbits), pattern):
        temp_word = vbits[i:i + pattern]
        temp_msg = 0
        for _ in range(pattern):
            try:
                temp_msg += temp_word[_] * (1 << _)
            except IndexError:
                # print('error:%d'%_)
                # print(len(vbits))
                pass
        msg.append(temp_msg)
    return msg

def wasserstein_dist(refer: List[float], input: List[float]):
    """
    计算两个正态分布之间的距离
    :param refer: 正态分布，即ATCG的比例关系
    :param input: 同ratio1
    :return: wasserstein_distance
    """
    p = np.array(refer)
    q = np.array(input)
    emd = wasserstein_distance(p, q)
    return emd

def js_dist(refer: List[float], input: List[float]):
    res = distance.jensenshannon(refer, input)
    return res

def calgc(ratio: dict):
    """
    计算GC含量
    :param ratio: ATCG比例
    :return: GC含量
    """
    gccontent = (ratio['G'] + ratio['C']) / sum(
        [ratio[l] for l in ratio.keys()])
    return gccontent

def calindex(ratio, resolution: int):
    """
    返回ratio对应的哈希值，用于统计是否访问过该ratio
    :param ratio: ATCG比例
    :param resolution: ATCG加起来的值
    :return: ratio对应的哈希值
    """
    return ratio['A'] * (resolution**3) + ratio['C'] * (
        resolution**2) + ratio['G'] * resolution + ratio['T']

class chn:

    def __init__(self, resolution, sigma, step):
        self.pattern = 8  # 8个二进制位对应一个letter
        self.step = step  # 每次滑窗步进的二进制位数，--> code rate
        self.pattern_mask = (1 << 7 * self.pattern) - 1  # 对应pattern的mask
        # self.pattern_mask = (1 << self.pattern) - 1  # 对应pattern的mask
        self.MAX_SEQ = 250000  # 允许的最大Hypothesis数量 250000
        self.PENALTY_THRESHOLD = 10  # 10 控制解码退出，惩罚太高直接退出
        self.NODE_THRESHOLD = 50 # 50  # 解码树每层最多节点数
        self.BUCKET_NUM = 13 # 10

        self.mod = 0
        self.duplicate = 1
        self.resolution = resolution
        self.sigma = sigma
        self.alphabet = []
        self.alphabet_set(self.resolution, self.sigma)

    def alphabet_set(self, resolution: int, sigma: int):
        """
        计算合法的全部letter, sigma mod
        ACGTMKRY
        :param resolution: ATCG加起来的值
        :return: None
        """
        base_A = composite_letter(code_pattern(resolution)['A'])
        base_C = composite_letter(code_pattern(resolution)['C'])
        base_G = composite_letter(code_pattern(resolution)['G'])
        base_T = composite_letter(code_pattern(resolution)['T'])
        base_M = composite_letter(code_pattern(resolution)['M'])
        base_K = composite_letter(code_pattern(resolution)['K'])
        base_R = composite_letter(code_pattern(resolution)['R'])
        base_Y = composite_letter(code_pattern(resolution)['Y'])

        if sigma == 4:
            self.alphabet = [base_A, base_C, base_G, base_T]
        elif sigma == 6:
            self.alphabet = [base_A, base_C, base_G, base_T, base_M, base_K]
        elif sigma == 8:
            self.alphabet = [
                base_A, base_C, base_G, base_T, base_M, base_K, base_R, base_Y
            ]
        self.mod = sigma
        

    def alphabet_show(self, n=10):
        """
        show alphabet
        :param n: 每行显示的letter个数
        :return: None
        """
        count = 0
        for i in range(self.mod):
            print(self.alphabet[i], end=' ')
            count += 1
            if count == n:
                print()
                count = 0
        print()
        print(f'resolution: {self.resolution}\t\tsigma: {self.sigma}')

    def dnacallowed(self, codetext):
        """
        Limit homopolymer length and GC content
        """
        all_base = ['A','C','G','T','M','K','R','Y']
        related_alphabet = [[[2,0,0,0],[1,1,0,0],[1,0,1,0]], # A,M,R
                            [[0,2,0,0],[1,1,0,0],[0,1,0,1]], # C,M,Y
                            [[0,0,2,0],[0,0,1,1],[1,0,1,0]], # G,K,R
                            [[0,0,0,2],[0,0,1,1],[0,1,0,1]]] # T,K,Y
        maxgcratio = 0.55 # 0.55
        mingcratio = 0.45 # 0.45
        max_runs = 4 # 4
        max_window = 160
        is_homopolymer = False
 
        alphabetletter = [[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 2],
                          [1, 1, 0, 0], [0, 0, 1, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
               
        alphabet_to_delete = []
        if len(codetext) > max_runs-1:
            for sub_related_alphabet in related_alphabet:
                tmp = []
                for item in codetext[-max_runs:]:
                    apb = [item.ratio_dict[x] for x in all_base[:4]]
                    tmp.append(apb in sub_related_alphabet)
                is_homopolymer = all(tmp)
                if is_homopolymer:
                    alphabet_to_delete = [alphabetletter.index(i) for i in sub_related_alphabet]

        window = len(codetext) if len(codetext)<= max_window else max_window
        if len(codetext) > max_runs-1:
            alphabet_count_in_window = [0] * 8
            for code_idx in range(-window , 0):
                current_code_base_ratio = [codetext[code_idx].ratio_dict[x] for x in all_base[:4]]
                for _ in range(len(alphabetletter)):
                    alphabet_count_in_window[_] += 1 if current_code_base_ratio == alphabetletter[_] else 0
            
            gc_count = (alphabet_count_in_window[1]+alphabet_count_in_window[2])
            my_count = (alphabet_count_in_window[4]+alphabet_count_in_window[7])
            kr_count = (alphabet_count_in_window[5]+alphabet_count_in_window[6])

            if mingcratio *len(codetext) > gc_count and self.sigma == 4:
                alphabet_to_delete += [0,3]
            if maxgcratio *len(codetext) < gc_count and self.sigma == 4:
                alphabet_to_delete += [1,2]
            
            if my_count < kr_count:
                alphabet_to_delete += [5,6]
            elif my_count > kr_count:
                alphabet_to_delete += [4,7]
            
            alphabet = self.alphabet
            alphabet = [item for index, item in enumerate(alphabet) if index not in list(set(alphabet_to_delete))]
            
            mod = len(alphabet)
        else:
            alphabet = self.alphabet[:self.sigma]
            mod = len(alphabet)
        
        return mod, alphabet


    def encode_c(self, vbits):
        """
        编码
        if DEBUG0: msg为全0数据
        if DEBUG1: msg为全1数据
        :param msg: input 数据
        :return: composite hedges编码后letter
        """
        flag_high = True
        if DEBUG0:
            vbits = [0] * 256
        if DEBUG1:
            vbits = [1] * 256
        nm = len(vbits)
        if nm > self.MAX_SEQ:
            raise ValueError('encode: MAXSEQ too small')
        codetext = []
        rawtext = []
        prev_bits = 0
        k = 0
        while k < nm:
            bit = 0
            for _ in range(self.step):
                try:
                    bit = (bit << 1) + vbits[k]
                    k += 1
                except:
                    bit = bit << 1
                    k += 1
            prev_bits = ((prev_bits << self.step) + bit) & self.pattern_mask

            mod, alphabet = self.dnacallowed(codetext)

            digest_code = self.digest(prev_bits, k - 1, mod)
            rawtext.append(digest_code)
            letter = alphabet[digest_code]
            temp = deepcopy(letter)
            codetext.append(temp)


        return codetext

    def decode_c(self, data: List[composite_letter]):
        """
        通过遍历每一步的全部可能性，计算可能生成的当前letter可能性与data比较，并返回概率最大的数据
        :param data: letter list
        :return: decode后的二进制数据及其的得分
        """
        hypotree = HypothesisTree(step=self.step)
        hypotree_plist = [hypotree.root]
        hypotree_qlist = []
        index = 0
        
        for _ in range(len(data)):
            hypotree_qlist = []
            index += self.step
            for hypotree_p in hypotree_plist:
                penalty_list = []
                hypo_prev_bits_list = []
                
                for bit in range(1 << self.step):                   
                    hypo_prev_bits = ((hypotree_p.prev_bits << self.step) +
                                      bit) & self.pattern_mask

                    mod, alphabet = self.dnacallowed(data[:_])
                    hypo_code = self.digest(hypo_prev_bits, index - 1,
                                            mod)
                    
                    hypo_letter = alphabet[hypo_code]

                    penalty = js_dist(hypo_letter.ratio_cal(),
                                      data[_].ratio_cal())
                    
                    penalty_list.append(penalty)
                    hypo_prev_bits_list.append(hypo_prev_bits)

                hypotree.add_children(node=hypotree_p,
                                    penalty_list=penalty_list,
                                    prev_bits_list=hypo_prev_bits_list)
                # delete high penalty node
                for child in hypotree_p.children:
                    if child.penalty < self.PENALTY_THRESHOLD:
                        hypotree_qlist.append(child)

            if len(hypotree_qlist) > self.NODE_THRESHOLD:
                bucket_list = [[] for _ in range(self.BUCKET_NUM)]
                bucket_size = self.PENALTY_THRESHOLD / self.BUCKET_NUM
                for i in range(len(hypotree_qlist)):
                    bucket_list[int(hypotree_qlist[i].penalty //
                                    bucket_size)].append(hypotree_qlist[i])
                cnt = 0
                templist = []
                min_bucket_size = len(bucket_list[0])
                for i in range(self.BUCKET_NUM):
                    templist.append(bucket_list[i])
                    cnt += len(bucket_list[i])
                    if cnt >= self.NODE_THRESHOLD:
                        break
                templist = [
                    element for sublist in templist for element in sublist
                ]
                hypotree_qlist = templist[:self.NODE_THRESHOLD
                                          if min_bucket_size < self.
                                          NODE_THRESHOLD else min_bucket_size]

            hypotree_plist = hypotree_qlist

        del hypotree_qlist
        return hypotree_plist

    def hypo_backtrack(self, hypotree_list):
        min_penalty = np.inf
        chosen = []
        res = []

        for hypo in hypotree_list:
            if hypo.penalty < min_penalty:
                min_penalty = hypo.penalty
        for i in range(len(hypotree_list)):
            if hypotree_list[i].penalty == min_penalty:
                chosen.append(i)
        for i in range(len(chosen)):
            stack = []
            stack_n = []
            p = hypotree_list[chosen[i]]
            while p.parent != None:
                stack.append(p.bits)
                p = p.parent
            stack.reverse()
            for num in stack:
                if self.step != 1:
                    stack_n.append((num // (self.step)) % self.step)
                    stack_n.append(num % self.step)
                else:
                    stack_n.append(num)
            res.append(stack_n)

        return res, min_penalty

    def digest(self, bits, k: int, mod: int):
        """
        计算并返回对应的哈希值
        :param bits: 前k位二进制数据
        :param k: 当前二进制位的index
        :param mod: 模
        :return:
        """
        return int(ran_hash(bits + (k << self.pattern)) % mod)

    def to_combine(self, encode_res: List[dict]) -> str:
        tempstr = ''
        cnt = 0
        for item in encode_res:
            for key in code_pattern(self.resolution):
                if item.ratio_dict == code_pattern(self.resolution)[key]:
                    tempstr += key
                    cnt += 1

        return tempstr

    def combine_to_letters(self, string: str) -> List[composite_letter]:
        data = []
        for i in range(len(string)):
            data.append(
                composite_letter(code_pattern(self.resolution)[string[i]]))
        return data