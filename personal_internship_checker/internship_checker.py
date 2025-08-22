import pandas as pd
import requests
from bs4 import BeautifulSoup
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
import schedule
import time

# Email Configuration
EMAIL_ADDRESS = 'andrew.patella04@gmail.com'
EMAIL_PASSWORD = 'kzxx ydkk iksd qcvs'


# Load URLs from Excel
def load_urls(file_path):
    df = pd.read_excel(file_path)
    return df[['Company', 'URL']].values.tolist()


# Check URLs for Keywords
def check_internships(urls, keywords):
    found_companies = []

    for company, url in urls:
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                soup = BeautifulSoup(response.text, 'html.parser')
                text = soup.get_text().lower()

                # Extract and process links/buttons correctly
                links = [a.get_text().strip().lower() for a in soup.find_all('a', href=True) if a.get_text()]
                buttons = [btn.get_text().strip().lower() for btn in soup.find_all(['button', 'input']) if
                           btn.get_text()]

                # Check keywords in links and buttons (correctly handling lists)
                for keyword in keywords:
                    if any(keyword.lower() in link for link in links) or \
                            any(keyword.lower() in button for button in buttons):
                        found_companies.append((company, url))
                        break  # Avoid adding duplicates if multiple matches are found
        except Exception as e:
            print(f"Error accessing {url}: {e}")

    return found_companies


# Send Email Notification
def send_email(companies):
    if companies:
        subject = 'Internship Opportunities Found'
        body = "The following companies have internship listings or mentions:\n\n"
        body += "\n".join([f" - {company} has internships available: {url}" for company, url in companies])
        body += "\n\n Happy hunting!"
    else:
        subject = 'No Internships Found'
        body = "No internship listings or mentions were found today."

    msg = MIMEMultipart()
    msg['From'] = EMAIL_ADDRESS
    msg['To'] = 'brendan.acosta@colorado.edu'
    msg['Subject'] = subject
    msg.attach(MIMEText(body, 'plain'))  # Send as plain text for simplicity

    with smtplib.SMTP('smtp.gmail.com', 587) as server:
        server.starttls()
        server.login(EMAIL_ADDRESS, EMAIL_PASSWORD)
        server.send_message(msg)


# Main Job to Run Daily
def job():
    urls = load_urls(r"C:/Users/ajpat/PycharmProjects/internshipChecker/company_urls.xlsx")
    keywords = ['intern', 'internship', 'summer intern', 'summer 2025']
    found_companies = check_internships(urls, keywords)
    send_email(found_companies)
    print("Internships compiled. Happy hunting!")


# Schedule to Run Daily at 9 AM
schedule.every().day.at("10:07").do(job)

# Run Immediately for Testing
job()

# Keep Script Running
while True:
    schedule.run_pending()
    time.sleep(60)