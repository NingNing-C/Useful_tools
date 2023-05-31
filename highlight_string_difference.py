import pandas as pd
import xlsxwriter

def highlight_chars(text, chars_to_highlight):
    highlighted_text = ''
    for char in text:
        if char in chars_to_highlight:
            highlighted_text += f'<font color="red">{char}</font>'
        else:
            highlighted_text += char
    return highlighted_text

# Example usage
data = {'String': ['Hello, how are you?', 'I love Python!', 'Data Science is awesome']}
df = pd.DataFrame(data)

chars_to_highlight = ['o', 'e', 'y']  # Characters to highlight

df['Highlighted String'] = df['String'].apply(lambda x: highlight_chars(x, chars_to_highlight))

# Create an Excel writer using pandas
excel_writer = pd.ExcelWriter('highlighted_strings.xlsx', engine='xlsxwriter')
df.to_excel(excel_writer, sheet_name='Sheet1', index=False)

# Get the workbook and the worksheet
workbook = excel_writer.book
worksheet = excel_writer.sheets['Sheet1']

# Add a format for the colored text
colored_text_format = workbook.add_format({'font_color': 'red'})

# Apply the formatting to the highlighted cells
worksheet.set_column('C:C', None, colored_text_format)  # Assuming 'Highlighted String' column is 'C'

# Save the workbook and close the Excel writer
excel_writer.save()
excel_writer.close()
