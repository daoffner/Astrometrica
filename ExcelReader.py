import openpyxl
from openpyxl.cell import get_column_letter
wb = openpyxl.load_workbook('Possibles 3.xlsx', data_only=True)
sheet = wb.get_sheet_by_name('Sheet1')
x= get_column_letter(sheet.max_column)

#print (ord(x))       Prints the ascii value of the max column; needed for testing
x=str(chr(ord(x)))
y=str(sheet.max_row)
#print(y)               Prints the row size of the sheet; needed for testing
tuple (sheet['A1' : x+y])

for rowOfCellObjects in sheet['A1' : x+y]:
    for cellObj in rowOfCellObjects:
        print(cellObj.coordinate, cellObj.value)

    print("---------- END OF ROW --------------")

print(sheet['A2'].value)
