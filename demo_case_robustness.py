import os
import Chamaeleo
from Chamaeleo.CHN.CHNcodec import CHN
from Chamaeleo.utils.pipelines_mod import RobustnessPipeline


if __name__ == "__main__":
    root_path = os.path.dirname(Chamaeleo.__file__)

    file_paths = {
        # "sme logo.jpg": os.path.join(root_path, "data", "CHN_test", "sme logo.jpg")
        "sme introduction.txt": os.path.join(root_path, "data", "CHN_test", "sme introduction.txt")
    }

    coding_schemes = {
        # "Base": BaseCodingAlgorithm(), "Church et al.": Church()
        "CHN": CHN()
    }
    error_corrections = {
        # "None": None, "Hamming": Hamming(), "ReedSolomon": ReedSolomon()
        "None": None
    }

    needed_indices = [
        # True, True
        False
    ]

    pipeline = RobustnessPipeline(
        coding_schemes=coding_schemes,
        error_corrections=error_corrections,
        needed_indices=needed_indices,
        file_paths=file_paths,
        nucleotide_insertion=0.001,
        nucleotide_mutation=0.03,
        nucleotide_deletion=0.001,
        sequence_loss=0,
        iterations=1,
        segment_length=288,
        index_length=16,
        need_logs=True,
        rb = True
    )

    pipeline.evaluate()
    pipeline.output_records(type="string")
