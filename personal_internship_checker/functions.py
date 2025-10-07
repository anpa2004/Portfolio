import requests
from bs4 import BeautifulSoup
import pandas as pd
import pickle
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from urllib.parse import urljoin
from datetime import datetime
import os
from urllib.parse import urljoin, urlparse
import re
import time
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.common.keys import Keys
from webdriver_manager.chrome import ChromeDriverManager

# ===============================
# Step 1: Load Company List
# ===============================
def load_companies(csv_file):
    return pd.read_csv(csv_file)

# ===============================
# Step 2: Job Parsers
# ===============================

def parse_greenhouse(url):
    jobs = []
    r = requests.get(url, timeout=10)
    soup = BeautifulSoup(r.text, "html.parser")
    postings = soup.select("div.opening")  # Greenhouse structure
    for post in postings:
        title = post.find("a").text.strip()
        link = urljoin(url, post.find("a")["href"])
        location = post.find("span", class_="location")
        jobs.append({
            "title": title,
            "link": link,
            "location": location.text.strip() if location else None,
            "date": None,
            "salary": None
        })
    return jobs

def parse_lever(url):
    jobs = []
    r = requests.get(url, timeout=10)
    soup = BeautifulSoup(r.text, "html.parser")
    postings = soup.select("div.posting")
    for post in postings:
        title = post.find("h5").text.strip()
        link = post.find("a")["href"]
        location = post.find("span", class_="sort-by-location")
        jobs.append({
            "title": title,
            "link": link,
            "location": location.text.strip() if location else None,
            "date": None,
            "salary": None
        })
    return jobs


def parse_workday(url, max_jobs=200, wait=5, scroll_pause=3):
 
    options = Options()
    options.add_argument("--headless=new")   # remove if you want to see the browser
    options.add_argument("--disable-gpu")
    options.add_argument("--no-sandbox")

    driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=options)
    driver.get(url)
    time.sleep(wait)

    jobs = []
    seen = set()
    body = driver.find_element(By.TAG_NAME, "body")

    while len(jobs) < max_jobs:
        job_cards = driver.find_elements(By.CSS_SELECTOR, "li.css-1q2dra3")
        if not job_cards:
            job_cards = driver.find_elements(By.CSS_SELECTOR, "li[data-automation-id='compositeTile'], div[data-automation-id='jobPosting']")

        for card in job_cards:
            try:
                # Title: prefer <a> text, fallback to subtitle
                title = "N/A"
                link = url
                try:
                    a_el = card.find_element(By.CSS_SELECTOR, "a[href]")
                    title = a_el.text.strip() if a_el.text.strip() else title
                    link = a_el.get_attribute("href")
                except:
                    pass

                if title == "N/A":
                    try:
                        title_el = card.find_element(By.CSS_SELECTOR, "ul[data-automation-id='subtitle'] li")
                        title = title_el.text.strip()
                    except:
                        pass

                # Location
                try:
                    location_el = card.find_element(By.CSS_SELECTOR, "div[data-automation-id='locations']")
                    location = location_el.text.strip()
                except:
                    location = "N/A"

                job_key = (title, location, link)
                if job_key not in seen and title != "N/A":
                    jobs.append({
                        "title": title,
                        "location": location,
                        "link": link
                    })
                    seen.add(job_key)
            except Exception as e:
                print(f"[Workday] Could not parse a job card: {e}")
                continue

        # Scroll and wait for lazy load
        body.send_keys(Keys.END)
        time.sleep(scroll_pause)

        if len(jobs) >= max_jobs or len(job_cards) == len(seen):
            break

    driver.quit()
    return jobs




# ===============================
# Step 3: Parser Router
# ===============================

# REPLACE WITH SCRAPER #

# ===============================
# Step 4: Deduplicate with Pickle
# ===============================
def load_seen(file="seen.pkl"):
    if os.path.exists(file):
        with open(file, "rb") as f:
            return pickle.load(f)
    return set()

def save_seen(seen, file="seen.pkl"):
    with open(file, "wb") as f:
        pickle.dump(seen, f)

# ===============================
# Step 5: Keyword Filter
# ===============================
def filter_jobs(jobs, keywords):
    results = []
    for job in jobs:
        if any(kw.lower() in job["title"].lower() for kw in keywords):
            results.append(job)
    return results

# ===============================
# Step 6: Email Report
# ===============================
def send_email(new_jobs, recipient, sender, password):
    msg = MIMEMultipart("alternative")
    msg["Subject"] = "Daily Job Report"
    msg["From"] = sender
    msg["To"] = recipient

    html = "<h3>New Jobs Found</h3><ul>"
    for comp, jobs in new_jobs.items():
        html += f"<li><b>{comp}</b><ul>"
        for job in jobs:
            html += f'<li><a href="{job["link"]}">{job["title"]}</a> - {job.get("location","")}</li>'
        html += "</ul></li>"
    html += "</ul>"

    msg.attach(MIMEText(html, "html"))

    with smtplib.SMTP_SSL("smtp.gmail.com", 465) as server:
        server.login(sender, password)
        server.sendmail(sender, recipient, msg.as_string())

# ===============================
# Step 7: Main Runner
# ===============================


### Running it

