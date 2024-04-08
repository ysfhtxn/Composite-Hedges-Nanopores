import os
import CHN
from CHN.CHNcodec import chn
from CHN.pipelines_mod import RobustnessPipeline


if __name__ == "__main__":
    root_path = os.path.dirname(CHN.__file__)

    file_paths = {
        "sme introduction.txt": os.path.join(root_path, "CHN_test", "sme introduction.txt")
    }

    coding_schemes = {
        "CHN": chn()
    }
    error_corrections = {
        "None": None
    }

    needed_indices = [
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
