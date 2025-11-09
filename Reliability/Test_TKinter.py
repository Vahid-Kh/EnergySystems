import tkinter as tk
from tkinter import ttk
from tkinter import filedialog as fd
from tkinter.messagebox import showinfo
import pandas as pd

# create the root window
root = tk.Tk()
root.title('Tkinter Open File Dialog')
root.resizable(False, False)
root.geometry('300x150')


lca_data = 0

df_lca = pd.read_excel("https://danfoss.sharepoint.com/:x:/r/sites/RACTechRoadmap/Shared%20Documents/General/Sustainability%20Tool/220203%20-%20TEMPLATE%20-%20Design%20for%20Sustainability.xlsx?d=w67d49b02dd014e76b982b1e37dabe7a3&csf=1&web=1&e=F2F1Qd"
, sheet_name="A1 - Raw material", index_col=0)
#
def lca_calc(df):
    df.head()
    return


def lca_calc(df):
    df.head()
    return

def select_file():
    filetypes = (
        ('Excel file', '*.xlsx'),
        ('All files', '*.*')
    )
    filename = fd.askopenfilename(
        title='Open a file',
        initialdir='/',
        filetypes=filetypes)
    showinfo(
        title='Selected File',
        message=filename
    )
    df_lca = pd.read_excel(filename, sheet_name="A1 - Raw material", index_col=0)

    return filename


# open button
open_button = ttk.Button(
    root,
    text='Open Database LCA Excel File',
    command=select_file
)


open_button2 = ttk.Button(
    root,
    text='Open Should Costing File',
    command=select_file
)



open_button.pack(expand=True)
open_button2.pack(expand=True)

# run the application
root.mainloop()