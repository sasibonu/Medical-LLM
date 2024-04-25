from app import db
from sqlalchemy import Column, Integer, String, Time, ForeignKey, JSON
from datetime import time

class User(db.Model):
    __tablename__ = 'users'
    id = db.Column(Integer, primary_key=True)
    email = db.Column(String(120), unique=True, nullable=False)
    password = db.Column(String(162), nullable=False)

    # Relationships with other entities
    search_histories = db.relationship('SearchHistory', backref='user', lazy=True, cascade="all, delete-orphan")
    sentiment_entries = db.relationship('Sentiment', backref='user', lazy=True, cascade="all, delete-orphan")

    # Notification tokens
    ios_token = db.Column(String(256), nullable=True)
    android_token = db.Column(String(256), nullable=True)
    web_token = db.Column(String(256), nullable=True)
    windows_token = db.Column(String(256), nullable=True)
    mac_token = db.Column(String(256), nullable=True)
    linux_token = db.Column(String(256), nullable=True)
    browser_token = db.Column(String(256), nullable=True)
    email_notification_token = db.Column(String(256), nullable=True)

    # Preferred notification time and timezone
    notification_time = db.Column(Time, nullable=True)
    timezone = db.Column(String(255), nullable=True)
    
    stripetoken = db.Column(String(255), nullable=True)
    
class Sentiment(db.Model):
    __tablename__ = 'sentiment'
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('users.id'), nullable=False)
    timestamp = db.Column(db.DateTime, default=db.func.current_timestamp())
    tag_name = db.Column(String(20), nullable=True)
    confidence = db.Column(db.Double, nullable=True)
    

class SearchHistory(db.Model):
    __tablename__ = 'search_histories'
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('users.id'), nullable=False)
    search_query = db.Column(db.String(255), nullable=False)
    timestamp = db.Column(db.DateTime, default=db.func.current_timestamp())

    def __repr__(self):
        return f"<SearchHistory id={self.id} user_id={self.user_id} search_query={self.search_query}>"