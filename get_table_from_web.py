import requests
from bs4 import BeautifulSoup
import pandas as pd

# URL of the website with the table
# url = 'https://en.wikipedia.org/wiki/List_of_countries_by_population_(United_Nations)'
url ='https://www.worldometers.info/geography/how-many-countries-in-europe/'

# Send a GET request to the website and parse the HTML content using BeautifulSoup
response = requests.get(url)
soup = BeautifulSoup(response.text, 'html.parser')

# Find the table in the HTML content
table = soup.find('table')

# Extract the table rows and columns into a Pandas dataframe
rows = []
for row in table.find_all('tr'):
    rows.append([cell.get_text(strip=True) for cell in row.find_all(['th', 'td'])])
df = pd.DataFrame(rows[1:], columns=rows[0])

# Print the dataframe
print(df)