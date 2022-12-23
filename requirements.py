import sys
import subprocess
import pkg_resources
from pkg_resources import DistributionNotFound, VersionConflict


def should_install_requirement(requirement):
    should_install = False
    try:
        pkg_resources.require(requirement)
    except (DistributionNotFound, VersionConflict):
        should_install = True
    return should_install


def install_packages(requirement_list):
    try:
        requirements = [requirement for requirement in requirement_list]
        print(requirements)
        if len(requirements) > 0:
            subprocess.check_call(
                [sys.executable, "-m", "pip", "install", *requirements]
            )
        else:
            print("Requirements already satisfied.")

    except Exception as e:
        print(e)


requirement_list = [
    "sqlalchemy",
    "pandas",
    "pg8000==1.16.6",
    "tqdm",
    "PySide6",
    "matplotlib",
    "requests",
    "plotly",
    "Biopython",
]

print("Checking requirements.")
install_packages(requirement_list)
print("Done.")
