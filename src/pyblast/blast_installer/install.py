import platform
import subprocess
import requests
from urllib.parse import urljoin


BLAST_VERSION = "2.15.0"
BASE_URL = f"https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/"


def install_blast(blast_version: str):
    """
    Installs the latest version of BLAST based on the OS.

    Args:
    base_url (str): The base URL where the BLAST executables are hosted.

    Returns:
    str: Installation status message.
    """
    base_url = BASE_URL + blast_version + '/'
    os_system = platform.system()
    os_architecture = platform.machine()

    # Based on the OS and architecture, define the file to download
    if os_system == "Windows" and os_architecture == "AMD64":
        filename = f"ncbi-blast-{blast_version}+-x64-win64.tar.gz"
    elif os_system == "Linux":
        filename = f"ncbi-blast-{blast_version}+-x64-linux.tar.gz" if os_architecture == "x86_64" else f"ncbi-blast-{blast_version}+-aarch64-linux.tar.gz"
    elif os_system == "Darwin":
        filename = f"ncbi-blast-{blast_version}+-x64-macosx.tar.gz"
    else:
        return "Unsupported OS or architecture."

    # Construct the full URL for the file to download
    download_url = urljoin(base_url, filename)

    # Download the file
    response = requests.get(download_url)
    with open(filename, 'wb') as f:
        f.write(response.content)

    # For simplicity, we will assume Linux/Unix-like commands for extraction
    # Windows may require a different method or an additional dependency to handle tar.gz files
    subprocess.run(['tar', '-xzf', filename])
    return f"BLAST installed successfully from {filename}"

#install_status = install_blast(BLAST_VERSION)
#print(install_status)


def test_blast_installation(blast_executable_path):
    """
    Tests whether BLAST is installed by checking the executable path.
    
    Args:
    blast_executable_path (str): The path to the BLAST executable.
    
    Returns:
    bool: True if BLAST is correctly installed, False otherwise.
    """
    try:
        # Try to run the 'blastn -version' command and check if it succeeds
        result = subprocess.run([blast_executable_path, '-version'], capture_output=True, text=True, check=True)
        return 'blastn:' in result.stdout
    except subprocess.CalledProcessError:
        return False


