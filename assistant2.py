import nltk
from nltk.corpus import stopwords
from nltk.tokenize import word_tokenize

import requests
from bs4 import BeautifulSoup
from openai import OpenAI
from Bio import Entrez
import os

# Download NLTK resources (run this once)
# nltk.download('punkt')
# nltk.download('stopwords')

Entrez.email = "dhba507@colorado.edu"
def is_therapy_needed(user_input):
    stop_words = set(stopwords.words('english'))
    
    therapy_keywords = [
        "anxiety", "depression", "stress", "mental health", "therapy", "counseling", 
        "psychology", "psychiatry", "well-being", "emotional health", "self-care"
    ]
    
    therapy_keywords += list(stop_words)
    words = word_tokenize(user_input.lower())
    
    for word in words:
        if word in therapy_keywords:
            return True
    return False

def therapist_assistant(user_input):
    if is_therapy_needed(user_input):
        core_issue, therapy_aspect = identify_issue(user_input)
        print("YOUR ISSUE: ", core_issue)
        therapy_context = get_therapy_context(core_issue, therapy_aspect)
        conduct_therapy(therapy_context)
    else:
        print("It seems like therapy is not needed at the moment.")


# This function needs to be implemented. This is just a placeholder right now.
def identify_issue(user_input):
    core_issue = "anxiety" 
    therapy_aspect = "Cognitive Behavioral Therapy"  
    return core_issue, therapy_aspect

def get_therapy_context(core_issue, therapy_aspect):
    try:
        handle = Entrez.esearch(db="pubmed", term=core_issue, retmax=3)
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

        titles = [entry['title'] for entry in articles]
        therapy_context = {
        "message": "These are the relevant articles from NLM",
        "articles": {'\n\n'.join(titles)},
        "recommendations": ["Give some suggestion"],
        "context": "Actually what happened."
    }
        
    except Exception as e:
        print("Error occurred during NLM search:", e)
        return []
    

    return therapy_context

def conduct_therapy(therapy_context):
    print("Therapist: Based on your concern, here are some recommendations:")
    for recommendation in therapy_context["recommendations"]:
        print("- " + recommendation)
    print("Therapist: " + therapy_context["context"])
    print("articles: " ,therapy_context["articles"])

while True:
    user_input = input("You: ")
    if user_input.lower() == "exit":
        print("Exiting therapy session.")
        break
    therapist_assistant(user_input)
