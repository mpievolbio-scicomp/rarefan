from flask import render_template
from app.views import LoginForm
from app import app

@app.route('/')
def index():
    return render_template("index.html")


@app.route('/login')
def login():
    form = LoginForm()
    return render_template('login.html', title='Sign In', form=form)
