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
import json

app = Flask(__name__)

# Download NLTK resources (run this once)
# nltk.download('punkt')
# nltk.download('stopwords')

Entrez.email = "dhba507@colorado.edu"

# Connect to Redis
redis_client = redis.StrictRedis(host='localhost', port=6379, db=0)

@app.route("/therapy", methods=["POST"])
def therapist_endpoint():
    print("Starting to execute")
    data = request.get_json()
    user_input = data["user_input"]
    core_issue, therapy = is_therapy_needed(user_input, mesh)
    print("Detected Issue: ", core_issue)
    if therapy:
        therapy_context = get_therapy_context(core_issue)
        result = {
            "Your Issue": core_issue, 
            "Is Therapy Needed": True, 
            "Therapy Context": therapy_context
        }
        # return jsonify({"Your Issue": core_issue, "Is Therapy Needed": True, "Therapy Context": therapy_context})
        return result
    else:
        return jsonify({"Is Therapy Needed": "It seems therapy is not needed at the moment!"})

def is_therapy_needed(user_input, mesh_data):
    names = []
    for key, value in mesh_data.items():
        if 'name' in value:
            names.append(value['name'])

    stop_words = set(stopwords.words('english'))

    therapy_keywords = names

    therapy_keywords = [word for word in therapy_keywords if word not in stop_words]
    new_keyword = [word.lower() for word in therapy_keywords]
    words = word_tokenize(user_input.lower())
    for word in words:
        if word in new_keyword:
            return word, True
    return None, False

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
            article_info['title'] = summary_result[0]['Title']
            article_info['url'] = f"https://pubmed.ncbi.nlm.nih.gov/{article_id}/"
            article_info['pmid'] = article_id

            if summary:
                get_result = summary.decode("utf-8")
                article_info['summary'] = get_result
            else:
                article_info['summary'] = scrape_content(article_info['url'])
                redis_client.set(article_id, article_info['summary'])
            end_time = time.time()

            print(f"Time taken to fetch summary for article {article_id}: {end_time - start_time} seconds")
            articles.append(article_info)
            

        therapy_context = []
        for article_info in articles:
            context = {
                "Suggestions": {
                    "Article Title": article_info['title'],
                    "Article pmid": article_info['pmid'],
                    "Article Link": article_info['url'],
                    "Article Summary": article_info['summary']
                }
            }
            therapy_context.append(context)       

    except Exception as e:
        print("Error occurred during NLM search:", e)
        return {"recommendations": ["No recommendations available due to an error"]}

    return therapy_context

def scrape_content(url):
    print("Scraping started")
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')
    content = ' '.join([p.text for p in soup.find_all('p')])
    return generate_summary(content)

def generate_summary(content):
    tokenizer = BartTokenizer.from_pretrained('facebook/bart-large-cnn')
    model = BartForConditionalGeneration.from_pretrained('facebook/bart-large-cnn')
    print("Summary BART running")
    inputs = tokenizer.encode("summarize: " + content, return_tensors="pt", max_length=1024, truncation=True)
    summary_ids = model.generate(inputs, num_beams=4, max_length=150, early_stopping=True)
    summary = tokenizer.decode(summary_ids[0], skip_special_tokens=True)
    return summary

if __name__ == "__main__":
    app.run(debug=True)
