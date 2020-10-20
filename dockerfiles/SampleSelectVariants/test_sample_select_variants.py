import unittest
import subprocess
import os
import zarr


class TestScript(unittest.TestCase):

    def setUp(self):
        if os.path.exists("output.vcf"):
            os.remove("output.vcf")

    def test(self):
        result = subprocess.run(["python", "sample_select_variants.py",
                                 "--sample-genotypes", "fixture/sample_genotypes.zarr.zip",
                                 "--sites-called", "fixture/sites_called.zarr.zip",
                                 "--sites-selected", "fixture/sites_selected.zarr.zip",
                                 "--output", "output.vcf",
                                 "--contig", "2R",
                                 ],
                                check=True,
                                capture_output=True)
        print(result.stdout.decode())

        with open("output.vcf", mode="r") as output, \
             open("fixture/expected.vcf", mode="r") as expected:
            output_lines = output.readlines()
            expected_lines = expected.readlines()
            for el, ol in zip(expected_lines, output_lines):
                if el.startswith("##"):
                    # compare headers
                    key = el[2:].split("=")[0]
                    if key in {"fileDate", "sites_called", "sites_selected", "sample_genotypes"}:
                        # for these headers, just check the key is correct, value may differ
                        # depending on which system the test is run on
                        assert ol.startswith(f"##{key}="), (el, ol)
                    else:
                        assert el == ol, (el, ol)
                else:
                    assert el == ol, (el, ol)


if __name__ == '__main__':
    unittest.main()
