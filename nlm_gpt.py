
### README
### This code searches the NLM and returns the top 3 best-match articles for the given term. 
### The code also downloads the returned articles as pdf and txt files. They are successfully getting downloaded to my local.
### But the downloaded pdf files are not opening whereas the txt files are.
### Another function uses openAI to rank the articles but the number of requests is very low for free plan and I exhausted it. 


import requests
from bs4 import BeautifulSoup
from openai import OpenAI
from Bio import Entrez
import os

llm = OpenAI(api_key='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')
Entrez.email = "dhba507@colorado.edu"

def search_nlm(search_term):
    try:
        handle = Entrez.esearch(db="pubmed", term=search_term, retmax=3)
        search_results = Entrez.read(handle)
        handle.close()

        article_ids = search_results['IdList']
        articles = []

        for article_id in article_ids:
            article_info = {}
            summary_handle = Entrez.esummary(db="pubmed", id=article_id)
            summary_result = Entrez.read(summary_handle)
            summary_handle.close()

            article_info['title'] = summary_result[0]['Title']
            article_info['url'] = f"https://pubmed.ncbi.nlm.nih.gov/{article_id}/"
            article_info['pmid'] = article_id
            articles.append(article_info)

        return articles

    except Exception as e:
        print("Error occurred during NLM search:", e)
        return []


def rank_articles_with_gpt3(articles):
    try:
        prompt = "Rank the following articles based on relevance:\n\n"
        for idx, article in enumerate(articles):
            prompt += f"{idx+1}. {article['title']}\n{article['url']}\n\n"

        # response = llm.create_completion(
        #     model="gpt-3.5-turbo-instruct", 
        #     prompt=prompt,
        #     max_tokens=100,
        #     n=3,  # Request top 3 completions
        #     stop=None  # Do not stop generating completions
        # )
            
        response = llm.completions.create(
        model="gpt-3.5-turbo-instruct",
        prompt=prompt,
        max_tokens=30000,
        temperature=0.5,
        top_p=1,
        frequency_penalty=0,
        presence_penalty=0
    )

        ranked_articles = []
        for choice in response.choices:
            ranked_articles.append(articles[int(choice['text']) - 1])

        return ranked_articles
    
    except Exception as e:
        print("Error occurred during GPT-3 ranking:", e)
        return []
    
def download_articles_as_pdf(articles, folder_path):
    i=1
    for article in articles:
        try:
            pdf_url = f"https://pubmed.ncbi.nlm.nih.gov/{article['pmid']}/pdf/"
            response = requests.get(pdf_url)
            if response.status_code == 200:
                file_path = os.path.join(folder_path, f"{i}_{article['title']}.pdf")
                print("file-path:",file_path)
                with open(file_path, "wb") as f:
                    f.write(response.content)
                # print(f"Downloaded {article['title']} as PDF to {file_path}")
                print("Downloaded pdf file rank:",(i))
            else:
                print(f"Failed to download {article['title']}: PDF not available")
        except Exception as e:
            print(f"Error occurred while downloading {article['title']}: {e}")
        i=i+1


def download_articles_as_txt(articles, folder_path):
    i=1
    for article in articles:
        try:
            response = requests.get(article['url'])
            if response.status_code == 200:
                text_content = BeautifulSoup(response.text, 'html.parser').get_text()
                # text_content = ' '.join(text_content.split())
                text_content = '\n'.join(line.strip() for line in text_content.split('\n'))
                file_path = os.path.join(folder_path, f"{i}_{article['title']}.txt")
                with open(file_path, "w", encoding="utf-8") as f:
                    f.write(text_content)
                # print(f"Downloaded {article['title']} as TXT to {file_path}")
                print("Downloaded Text file rank:",(i))
            else:
                print(f"Failed to download {article['title']}: Article not available")
        except Exception as e:
            print(f"Error occurred while downloading {article['title']}: {e}")
        i=i+1



search_term = "diabetes"
articles = search_nlm(search_term)
print("Total number:",len(articles))
# for article in articles:
#     print(article,"\n")
folder_path = "/Users/dharinibaskaran/Desktop/Big Data/pdfs" # Edit path here
   


for article in articles:
    print(article,"\n")
download_articles_as_txt(articles, folder_path) 
download_articles_as_pdf(articles, folder_path)




# Commenting this because i reached my GPT request limits

# print("TOP 3 GPT Articles:")
# ranked_articles = rank_articles_with_gpt3(articles)
# for idx, article in enumerate(ranked_articles[:3]):
#     print(f"Ranked Article {idx+1}: {article['title']}")
#     print(f"URL: {article['url']}")
#     print()
