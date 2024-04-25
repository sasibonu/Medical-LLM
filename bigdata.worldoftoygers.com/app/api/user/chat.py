from flask import Blueprint, jsonify, request
# from flask_jwt_extended import jwt_required, get_jwt_identity
from config import Config
from app import db
import openai
import os
import time as tm
import json
from Bio import Entrez
import requests
from bs4 import BeautifulSoup

#########----------------------------------------------------------#############
# fetching details from config file

entrez_access_email_id = Config.entrez_access_email_id
#print("!!!!!!!entrez_access_email_id!!!!", entrez_access_email_id)
assistant_id_config = Config.assistant_id
#print("!!!!!!!assistant_id!!!!", assistant_id_config)
OPENAI_API_KEY_config = Config.OPENAI_API_KEY
#print("!!!!!!!OPENAI_API_KEY!!!!", OPENAI_API_KEY_config)

#########----------------------------------------------------------#############

Entrez.email = entrez_access_email_id
bp = Blueprint('chat', __name__)



def scrape_content(url):
    print("Scraping started")
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')
    content = ' '.join([p.text for p in soup.find_all('p')])
    return content

def get_psychology_context(issue):
    try:
        try:
            handle = Entrez.esearch(db="pubmed", term=issue, retmax=3)
            #print("!!!!!!!!!handle!!!!!!!!!!", handle)
        except Exception as e:
            print("Failed at esearch:", e)
            return [{'pmid': '1', 'link': '2', 'content': '3'}]
        try:
            search_results = Entrez.read(handle)
        except Exception as e:
            print("Failed at Entrez.read: ", e)
            return {"Error": "Failed at search results"}
        try:
            handle.close()
        except Exception as e:
            return {"Error": "Failed at handle close"}

        article_ids = search_results['IdList']
        print(article_ids)
        articles = []

        for article_id in article_ids:
            article_info={}
            summary_handle = Entrez.esummary(db="pubmed", id=article_id)
            summary_result = Entrez.read(summary_handle)
            summary_handle.close()

            article_info['pmid'] = article_id
            article_info['link'] = f"https://pubmed.ncbi.nlm.nih.gov/{article_id}/"
            article_info['content'] = scrape_content(article_info['link'])
            #print("!!!!!!!!before article!!!!!!!!!!!", article_info)
            articles.append(article_info)
            #print("!!!!!!!!after article!!!!!!!!!!!", article_info)
        
        #print("PMIDs:", articles)

    except Exception as e:
        print("Error occurred during NLM search:", e)
        return {"recommendations": ["No recommendations due to error"]}
    return articles

@bp.route('/chat', methods=['POST'])
# @jwt_required()
def interact_with_openai():
    # user_id = get_jwt_identity()  # Use JWT to authenticate the user
    data = request.get_json()
    user_message = data.get('message')
    if not user_message:
        return jsonify({"error": "No message provided"}), 400

    # OpenAI API interaction
    try:
        # Set your OpenAI API key here
        os.environ["OPENAI_API_KEY"] = OPENAI_API_KEY_config 
        
        # Hardcoded Assistant ID
        assistant_id = assistant_id_config

        client = openai.OpenAI()
        assistant = client.beta.assistants.retrieve(assistant_id)

        # If thread_id is not provided in request, create a new thread
        thread_id = data.get('thread_id', None)
        if not thread_id:
            thread = client.beta.threads.create()
            thread_id = thread.id

        # Create a message in the existing or new thread
        client.beta.threads.messages.create(
            thread_id=thread_id,
            role="user",
            content=user_message
        )
        
        # Create a run to process the message
        run = client.beta.threads.runs.create(
            thread_id=thread_id,
            assistant_id=assistant.id,
        )
        
        # Wait for the run to complete
        while True:
            tm.sleep(1)
            run = client.beta.threads.runs.retrieve(
                thread_id=thread_id,
                run_id=run.id
            )
            if run.status == "requires_action":
                query = (run.required_action.submit_tool_outputs.tool_calls)[0].function.arguments
                issue = json.loads(query)
                print(issue['query'])
                result = get_psychology_context(issue)
                #print("!!!!!!!!!Result!!!!!!!!!!", result)
                content_list = [entry['content'] for entry in result]
                page_content = ', '.join(content_list)
                run = client.beta.threads.runs.submit_tool_outputs(
  thread_id=run.thread_id,
  run_id=run.id,
  tool_outputs=[
    {
      "tool_call_id": run.required_action.submit_tool_outputs.tool_calls[0].id,
      "output": page_content
    }
  ]
)
            if run.status == "completed":
                break

        # Retrieve the latest message
        messages = client.beta.threads.messages.list(
            limit=1,
            order="desc",
            thread_id=thread_id
        )
        gpt_return = messages.data[0].content[0].text.value
        print("!!!!!!!!!!!!!!!!!!!", gpt_return)
    except Exception as e:
        return jsonify({"error": str(e)}), 500

    return jsonify({"response": gpt_return, "thread_id": thread_id})
