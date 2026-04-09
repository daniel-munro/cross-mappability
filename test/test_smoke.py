import os
import subprocess
import unittest


class WorkflowSmokeTest(unittest.TestCase):
    def test_toy_workflow_completes(self):
        env = os.environ.copy()
        env["TMPDIR"] = "/tmp"
        env["XDG_CACHE_HOME"] = "/tmp"
        required = [
            ["Rscript", "-e", "library(argparser); library(data.table); library(intervals); library(seqinr); library(stringr)"],
            ["bowtie", "--version"],
            ["bowtie-build", "--version"],
        ]
        for command in required:
            subprocess.run(
                command,
                env=env,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True,
            )


if __name__ == "__main__":
    unittest.main()
