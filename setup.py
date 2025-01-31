from setuptools import setup, find_packages
import os

readme_path = os.path.join(os.path.dirname(__file__), "README.md")
with open(readme_path, "r") as readme_file:
    readme = readme_file.read()

with open('requirements.txt', 'r') as f:
    required = f.read().splitlines()

with open("skeleton_shollstice/__init__.py") as f:
    for line in f:
        if line.startswith("__version__"):
            version = line.split('"')[1]

setup(
    name = 'skeleton_shollstice',
    version = version,
    description = "Shedding light on morphology skeletons through sholl analysis ",
    long_description=readme,
    author = "Matthew Mallory",
    author_email = "matt.mallory@alleninstitute.org",
    url = 'https://github.com/AllenInstitute/morpho_shollstice',
    packages = find_packages(),
    install_requires = required,
    include_package_data=True,
    # package_data={"morpho_shollstice": ["data/*"]},
    # setup_requires=['pytest-runner'],
)