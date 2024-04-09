import os
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from bs4 import BeautifulSoup
from urllib.parse import urljoin, urlparse

def requests_session(retries=3, backoff_factor=0.3, status_forcelist=(500, 502, 504)):
    session = requests.Session()
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session

def download_file(url, base_folder='downloads'):
    try:
        session = requests_session()  # Use the session with retries
        path_parts = urlparse(url).path.split('/')
        folder_structure = base_folder if len(path_parts) < 3 else os.path.join(base_folder, path_parts[-3], path_parts[-2])
        os.makedirs(folder_structure, exist_ok=True)
        response = session.get(url)  # Use session to make a request
        if response.status_code == 200:
            filename = os.path.join(folder_structure, os.path.basename(urlparse(url).path))
            with open(filename, 'wb') as f:
                f.write(response.content)
            print(f"Downloaded: {url}")
        else:
            print(f"Failed to download: {url}")
    except Exception as e:
        print(f"Error downloading {url}: {str(e)}")

def scrape_and_save(url):
    try:
        session = requests_session()
        response = session.get(url)
        if response.status_code == 200:
            soup = BeautifulSoup(response.text, 'html.parser')
            for link in soup.find_all('a', href=True):
                if link['href'].endswith('.pdf'):
                    absolute_url = urljoin(url, link['href'])
                    download_file(absolute_url)
    except Exception as e:
        print(f"Error scraping {url}: {str(e)}")

def iterate_directories(base_url):
    for i in range(100):
        dir1 = format(i, '02x')
        for j in range(100):
            dir2 = format(j, '02x')
            target_url = f"{base_url}/{dir1}/{dir2}/"
            print(f"Scraping: {target_url}")
            scrape_and_save(target_url)

if __name__ == "__main__":
    base_url = "https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_pdf"
    iterate_directories(base_url)
