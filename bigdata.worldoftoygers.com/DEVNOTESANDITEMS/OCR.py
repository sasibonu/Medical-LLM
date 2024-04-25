import sys
import fitz  # PyMuPDF

def extract_text_from_pdf(pdf_path):
    # Open the PDF file
    with fitz.open(pdf_path) as doc:
        text = ""
        # Iterate through each page
        for page in doc:
            # Extract text from the page and add it to the text variable
            text += page.get_text()

    # Determine the name for the output txt file
    txt_filename = pdf_path.rsplit('.', 1)[0] + '.txt'

    # Save the extracted text to a txt file
    with open(txt_filename, 'w', encoding='utf-8') as f:
        f.write(text)
    print(f"Text extracted and saved to {txt_filename}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py file.pdf")
    else:
        pdf_path = sys.argv[1]
        extract_text_from_pdf(pdf_path)
