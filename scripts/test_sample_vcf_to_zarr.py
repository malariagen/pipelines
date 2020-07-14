import unittest
import subprocess
import shutil
import os
import zarr


class TestScript(unittest.TestCase):

    def setUp(self):
        shutil.rmtree("output", ignore_errors=True)
        os.mkdir("output")

    def check_example_callset(self, callset):
        assert "NA00001" in callset
        for contig in "19", "20":
            assert contig in callset["NA00001"]
            assert "variants" in callset[f"NA00001/{contig}"]
            assert "DP" in callset[f"NA00001/{contig}/variants"]
            assert "calldata" in callset[f"NA00001/{contig}"]
            assert "GT" in callset[f"NA00001/{contig}/calldata"]
            assert "GQ" in callset[f"NA00001/{contig}/calldata"]
        assert (2,) == callset["NA00001/19/variants/DP"].shape
        assert (6,) == callset["NA00001/20/variants/DP"].shape
        assert (2, 1, 2) == callset["NA00001/19/calldata/GT"].shape
        assert (6, 1, 2) == callset["NA00001/20/calldata/GT"].shape
        assert (2, 1) == callset["NA00001/19/calldata/GQ"].shape
        assert (6, 1) == callset["NA00001/20/calldata/GQ"].shape

    def test_conversion(self):
        result = subprocess.run(["python", "sample_vcf_to_zarr.py",
                                 "--input", "fixture/example.vcf.gz",
                                 "--output", "output/example.zarr",
                                 "--sample", "NA00001",
                                 "--contig", "19",
                                 "--contig", "20",
                                 "--field", "variants/DP",
                                 "--field", "calldata/GT",
                                 "--field", "calldata/GQ",
                                 ],
                                check=True,
                                capture_output=True)
        print(result.stdout.decode())
        callset = zarr.open("output/example.zarr")
        self.check_example_callset(callset)

    def test_conversion_zip(self):
        result = subprocess.run(["python", "sample_vcf_to_zarr.py",
                                 "--input", "fixture/example.vcf.gz",
                                 "--output", "output/example.zarr",
                                 "--sample", "NA00001",
                                 "--contig", "19",
                                 "--contig", "20",
                                 "--field", "variants/DP",
                                 "--field", "calldata/GT",
                                 "--field", "calldata/GQ",
                                 "--zip"
                                 ],
                                check=True,
                                capture_output=True)
        print(result.stdout.decode())
        callset = zarr.open("output/example.zarr.zip")
        self.check_example_callset(callset)


if __name__ == '__main__':
    unittest.main()
