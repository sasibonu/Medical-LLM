import redis
import json
redis_client = redis.StrictRedis(host='localhost', port=6379, db=1)
link = []
linkid=1
redis_client.set(linkid, "summary1")
summary = redis_client.get(linkid)
print("About to fetch summary from Redis, if present")
print("summary: \n", summary)
# if summary is None:
#                 # Summary not found in cache, fetch from PubMed and scrape content
#     print("entering if")
#     article_info['title'] = summary_result[0]['Title']
#     article_info['url'] = f"https://pubmed.ncbi.nlm.nih.gov/{article_id}/"
#     article_info['pmid'] = article_id
#     article_info['summary'] = scrape_content(article_info['url'])
#                 # Store summary in cache
#     redis_client.set(article_id, article_info['summary'])
# else:
#                 # Summary found in cache, retrieve from Redis
#     print("entering else")
#     article_info['summary'] = summary.decode('utf-8')


from flask import Flask, request, jsonify
from nltk.corpus import stopwords
from nltk.tokenize import word_tokenize
# from rake_nltk import Rake
import requests
from bs4 import BeautifulSoup
from transformers import BartForConditionalGeneration, BartTokenizer
from Bio import Entrez
import redis
import time
from meshdata import mesh

app = Flask(__name__)

# Download NLTK resources (run this once)
# nltk.download('punkt')
# nltk.download('stopwords')

Entrez.email = "dhba507@colorado.edu"

def get_therapy_context(core_issue):
    try:
        handle = Entrez.esearch(db="pubmed", term=core_issue, retmax=3)
        search_results = Entrez.read(handle)
        handle.close()

        article_ids = search_results['IdList']
        articles = []

        for article_id in article_ids:
            article_info = {}
            start_time = time.time() 
            summary_handle = Entrez.esummary(db="pubmed", id=article_id)
            summary_result = Entrez.read(summary_handle)
            summary_handle.close()
            summary = redis_client.get(article_id)
            # print("Retrieved Details:")
            # print(get_result[0]['Title'])

            if summary:
                get_result = json.loads(summary)
                article_info['title'] = get_result[0]['Title']
                article_info['url'] = f"https://pubmed.ncbi.nlm.nih.gov/{article_id}/"
                article_info['pmid'] = article_id
            else:
                article_info['title'] = summary_result[0]['Title']
                article_info['url'] = f"https://pubmed.ncbi.nlm.nih.gov/{article_id}/"
                article_info['pmid'] = article_id
                # article_info['summary'] = scrape_content(article_info['url'])
                summary_json = json.dumps(summary_result)
                redis_client.set(article_id, summary_json)
            end_time = time.time()  # End measuring time
            print(f"Time taken to fetch summary for article {article_id}: {end_time - start_time} seconds")
            articles.append(article_info)

        therapy_context = []
        for article_info in articles:
            context = {
                "Suggestions": {
                    "Article Title": article_info['title'],
                    "Article pmid": article_info['pmid'],
                    "Article Link": article_info['url'],
                    # "Article Summary": article_info['summary']
                }
            }
            therapy_context.append(context)
        print("suggestions result:")
        print(therapy_context)

    except Exception as e:
        print("Error occurred during NLM search:", e)
        return {"recommendations": ["No recommendations available due to an error"]}
    return therapy_context


def main():
    print("Welcome to our chatbot! Type 'exit' when you want to end.")
    while True:
        user_input = input("You: ")
        if user_input.lower() == "exit":
            print("Exiting therapy session.")
            break
        get_therapy_context(user_input)

if __name__ == "__main__":
    main()



    # [{'Item': [], 
    #   'Id': '38643251', 
    #   'PubDate': '2024 Apr 20', 
    #   'EPubDate': '2024 Apr 20', 
    #   'Source': 'Sci Rep', 
    #   'AuthorList': ['Gurbuz E', 'Riby DM', 'South M', 'Hanley M'], 
    #   'LastAuthor': 'Hanley M', 
    #   'Title': 'Associations between autistic traits, depression, social anxiety and social rejection in autistic and non-autistic adults.', 
    #   'Volume': '14', 
    #   'Issue': '1', 
    #   'Pages': '9065', 
    #   'LangList': ['English'], 
    #   'NlmUniqueID': '101563288', 
    #   'ISSN': '', 
    #   'ESSN': '2045-2322', 
    #   'PubTypeList': ['Journal Article'],
    #     'RecordStatus': 'PubMed - in process', 
    #     'PubStatus': 'epublish', 
    #     'ArticleIds': {'pubmed': ['38643251'], 'medline': [], 'doi': '10.1038/s41598-024-59532-3', 'pii': '10.1038/s41598-024-59532-3'}, 
    #     'DOI': '10.1038/s41598-024-59532-3', 
    #     'History': {'pubmed': ['2024/04/21 00:42'], 'medline': ['2024/04/21 00:42'], 'received': '2023/08/17 00:00', 'accepted': '2024/04/11 00:00', 'entrez': '2024/04/20 23:26'}, 'References': [], 'HasAbstract': IntegerElement(1, attributes={}), 
    #     'PmcRefCount': IntegerElement(0, attributes={}), 
    #     'FullJournalName': 'Scientific reports', 
    #     'ELocationID': 'doi: 10.1038/s41598-024-59532-3', 
    #     'SO': '2024 Apr 20;14(1):9065'}]