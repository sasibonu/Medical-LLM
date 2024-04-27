from flask import Flask, request, jsonify
from textblob import TextBlob


app = Flask(__name__)

@app.route('/')
def home():
    return 'Hello, welcome to the chatbot!'

@app.route('/predict', methods=['POST'])
def predict():
    text = request.form['text']
    blob = TextBlob(text)
    sentiment = blob.sentiment.polarity
    if sentiment > 0.5:
        response = 'Positive'
    elif sentiment < -0.5:
        response = 'Negative'
    else:
        response = 'Neutral'
    return jsonify({'response': response})
