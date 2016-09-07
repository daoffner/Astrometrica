import openpyxl
wb = openpyxl.load_workbook('Possibles 3.xlsx', data_only=True)
sheet = wb.get_sheet_by_name('Sheet1')
tuple (sheet['A1' : 'R5661'])

for rowOfCellObjects in sheet['A1' : 'R5661']:
    for cellObj in rowOfCellObjects:
        print(cellObj.coordinate, cellObj.value)

    print("---------- END OF ROW --------------")

print(sheet['A2'].value)
