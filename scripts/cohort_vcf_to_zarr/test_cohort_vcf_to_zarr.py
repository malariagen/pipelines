import unittest
import subprocess
import shutil
import os
import zarr
import numpy as np


class TestScript(unittest.TestCase):

    def setUp(self):
        shutil.rmtree("output", ignore_errors=True)
        os.mkdir("output")

    def test_phased(self):

        # run conversion
        result = subprocess.run(["python", "cohort_vcf_to_zarr.py",
                                 "--input", "fixture/phased.vcf",
                                 "--output", "output/phased.zarr",
                                 "--contig", "2R",
                                 "--field", "variants/POS",
                                 "--field", "variants/REF:S1",
                                 "--field", "variants/ALT:S1",
                                 "--field", "variants/AC",
                                 "--field", "variants/AF",
                                 "--field", "variants/CM",
                                 "--field", "calldata/GT",
                                 "--alt-number", "1",
                                 ],
                                check=True,
                                capture_output=True)
        print(result.stdout.decode())

        # open and check output
        callset = zarr.open("output/phased.zarr")
        assert "2R" in callset
        contig_callset = callset["2R"]
        assert "variants" in contig_callset
        assert "calldata" in contig_callset
        variants = contig_callset["variants"]
        for f in "POS", "REF", "ALT", "AC", "AF", "CM":
            assert f in variants, f
            assert (190,) == variants[f].shape
            if f in {"POS", "AC"}:
                assert variants[f].dtype.kind == "i"
            if f in {"REF", "ALT"}:
                assert variants[f].dtype == "S1"
            if f in {"AF", "CM"}:
                assert variants[f].dtype.kind == "f"
        calldata = contig_callset["calldata"]
        assert "GT" in calldata
        gt = calldata["GT"][:]
        assert gt.shape == (190, 4, 2)
        assert gt.dtype == "i1"
        assert np.all(gt >= 0)
        assert np.all(gt <= 1)
        assert np.any(gt == 0)
        assert np.any(gt == 1)
        x = variants['POS'][:]
        assert [135, 168, 170] == variants['POS'][:3].tolist()
        assert [b'A', b'A', b'A'] == variants['REF'][:3].tolist()
        assert [b'T', b'G', b'T'] == variants['ALT'][:3].tolist()
        assert [0, 0, 0, 0, 1, 0, 0, 0] == gt[10].flatten().tolist()
        assert [0, 0, 0, 0, 0, 1, 0, 0] == gt[12].flatten().tolist()


if __name__ == '__main__':
    unittest.main()
