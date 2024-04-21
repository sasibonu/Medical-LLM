2 files are added to the "therapy" folder.

curl:
curl --location 'http://127.0.0.1:5000/therapy' \
--header 'Content-Type: application/json' \
--data '{
    "user_input": "I am having thoughts of depression"
}'


Sample output:
{
    "Is Therapy Needed": true,
    "Your Issue": "depression"
    "Therapy Context": [
        {
            "Suggestions": {
                "Article Link": "https://pubmed.ncbi.nlm.nih.gov/38643303/",
                "Article Summary": "The.gov means it’s official. Federal government websites often end in.gov or.mil. The https:// ensures that you are connecting to the official website and that any information you provide is encrypted and transmitted securely. Post-stroke arrhythmia might be an early distinguishable marker for the presence of PSD.",
                "Article Title": "Post-stroke arrhythmia could be a potential predictor for post-stroke depression.",
                "Article pmid": "38643303"
            }
        },
        {
            "Suggestions": {
                "Article Link": "https://pubmed.ncbi.nlm.nih.gov/38643251/",
                "Article Summary": "The.gov means it’s official. Federal government websites often end in.gov or.mil. The https:// ensures that you are connecting to the official website and that any information you provide is encrypted and transmitted securely. Autistic people frequently experience negative judgements from non-autistic people. Understanding responses to negative social judgement among autistic people is crucial.",
                "Article Title": "Associations between autistic traits, depression, social anxiety and social rejection in autistic and non-autistic adults.",
                "Article pmid": "38643251"
            }
        },
        {
            "Suggestions": {
                "Article Link": "https://pubmed.ncbi.nlm.nih.gov/38643242/",
                "Article Summary": "The.gov means it’s official. Federal government websites often end in.gov or.mil. The https:// ensures that you are connecting to the official website and that any information you provide is encrypted. The study shows that ICBT for alcohol misuse is associated with reduced drinking and comorbid mental health difficulties over time.",
                "Article Title": "Internet-delivered therapy for alcohol misuse: engagement, satisfaction, and outcomes when patients select their preference for therapist- or self-guided treatment.",
                "Article pmid": "38643242"
            }
        }
    ]
}


therapy_convo.ep:
This file has the endpoint for conducting therapy. 

redis_ep.py
This file has the same endpoint, I have added an extra redis cache for faster retrieval of the summary.

meshdata.py:
This file contains the list of medical terms called "mesh_terms", taken from the NLM library. This json is used to identify what issue the user is talking about.

import statement:
from transformers import BartForConditionalGeneration, BartTokenizer

I have used this import statement to run the existing BART model, which will read the article by its link and provide a short summary. This summary scraping takes 15-18 seconds per file. So, to reduce latency I have used the redis cache.