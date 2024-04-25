from flask import Blueprint, jsonify, request
import requests
from flask_jwt_extended import jwt_required, get_jwt_identity
from app import db
from app.models import SearchHistory

bp = Blueprint('search', __name__, url_prefix='/search')

@bp.route('', methods=['GET'])
@jwt_required()
def search_contracts():
    service_name = 'azure_search_service'
    index_name = 'index'
    api_key = 'key'  # Securely store this

    # Retrieve the search term from the query parameter
    search_term = request.args.get('term', '')
    filter_query = request.args.get('filter', '')

    # Construct the JSON payload with the user-provided search term
    data = {
        "count": True,
        "skip": 0,
        "top": 50,
        "searchMode": "any",
        "queryType": "semantic",
         "facets": [
            "State,count:5,sort:count",
            "ClassificationCode,count:5,sort:count",
            "City,count:5,sort:count",
            "Active,count:5,sort:count",
            "NaicsCode,count:5,sort:count",
            "Office,count:5,sort:count",
            "ZipCode,count:5,sort:count",
            "SetASideCode,count:5,sort:count",
            "CountryCode,count:5,sort:count",
            "SetASide,count:5,sort:count",
            "PostedDate,count:5,sort:count",
            "PopCountry,count:5,sort:count",
            "ResponseDeadLine,count:5,sort:count",
            "PopZip,count:5,sort:count",
            "PopState,count:5,sort:count",
            "Type,count:5,sort:count",
            "PopCity,count:5,sort:count",
            "BaseType,count:5,sort:count",
            "Department_Ind_Agency,count:5,sort:count",
            "AwardDate,count:5,sort:count",
            "ArchiveDate,count:5,sort:count",
            "ArchiveType,count:5,sort:count",
            "CGAC,count:5,sort:count",
            "OrganizationType,count:5,sort:count",
            "Awardee,count:5,sort:count"
        ],
        "search": search_term,
        "filter": filter_query
    }

    endpoint = f'https://{service_name}.search.windows.net/indexes/{index_name}/docs/search?api-version=2023-11-01'
    headers = {'Content-Type': 'application/json', 'api-key': api_key}

    response = requests.post(endpoint, headers=headers, json=data)
    
    current_user_id = get_jwt_identity()

    # Save the search query to the database
    new_search_history = SearchHistory(user_id=current_user_id, search_query=search_term)
    db.session.add(new_search_history)
    db.session.commit()
    
    return jsonify(response.json())

