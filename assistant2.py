import nltk
from nltk.corpus import stopwords
from nltk.tokenize import word_tokenize
from rake_nltk import Rake
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
        "psychology", "psychiatry", "well-being", "emotional health", "self-care", "asleep"
    ]
    
    therapy_keywords = [word for word in therapy_keywords if word not in stop_words]
    
    words = word_tokenize(user_input.lower())
    
    for word in words:
        if word in therapy_keywords:
            print("Inside the first def:", word)
            return word, True
    return None, False

def therapist_assistant(user_input):
    core_issue, therapy = is_therapy_needed(user_input)
    print(core_issue,therapy)
    if(therapy):
        therapy_context = get_therapy_context(core_issue)
        print(therapy_context)
        conduct_therapy(therapy_context)
    else:
        print("It seems like therapy is not needed at the moment.")
        
def identify_issue(user_input):

    rake = Rake()
    rake.extract_keywords_from_text(user_input)
    keywords = rake.get_ranked_phrases()
    print("COREISSUE:", keywords)
    therapy_aspect = "Cognitive Behavioral Therapy"  
    return keywords, therapy_aspect

def get_therapy_context(core_issue):    

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
        "articles": '\n\n'.join(titles),
        "recommendations": ["Give some suggestion"],
        "context": "Actually what happened."
    }
        
    except Exception as e:
        print("Error occurred during NLM search:", e)
        return {"recommendations": ["No recommendations due to error"]}
    

    return therapy_context

def conduct_therapy(therapy_context):
    print("Therapist: Based on your concern, here are some recommendations:")
    if "recommendations" in therapy_context:
        for recommendation in therapy_context["recommendations"]:
            print("- " + recommendation)
    else:
        print("No recommendations available.")
        
    print("Therapist: " + therapy_context["context"])
    print("Articles: " ,therapy_context["articles"])


def main():
    print("Welcome to our chatbot! Type 'exit' when you want to end.")
    while True:
        user_input = input("You: ")
        if user_input.lower() == "exit":
            print("Exiting therapy session.")
            break
        therapist_assistant(user_input)

if __name__ == "__main__":
    main()

