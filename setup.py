from setuptools import setup

setup (
  name="PINE",
  version="0.1",
  description="Protein Interaction Network Extractor",
  license="Apache License 2.0",
  author_email="GroupHeartBioinformaticsSupport@cshs.org",
  install_requires = [
    "certifi==2019.9.11",
    "chardet==3.0.4",
    "decorator==4.4.1",
    "idna==2.8",
    "networkx==2.4",
    "numpy==1.17.3",
    "pandas==0.25.2",
    "py2cytoscape==0.7.1",
    "pydot==1.4.1",
    "pydotplus==2.0.2",
    "pyparsing==2.4.2",
    "python-dateutil==2.8.0",
    "python-igraph==0.7.1.post6",
    "pytz==2019.3",
    "requests==2.22.0",
    "six==1.12.0",
    "urllib3==1.25.6",
  ],
  dependency_links=["https://github.com/fivebillionmph/py2cytoscape/tarball/pine#egg=py2cytoscape-0.7.1"],
)