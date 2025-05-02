from tkinter import *
from tkinter import ttk
import pandas as pd
""" ------------Create Drop down menu and when selected text appears -------- """
# Create object
# root = Tk()

# Adjust size
# root.geometry("200x200")


# Change the label text
# def show():
#     label.config(text=clicked.get())

# # """------------------------- " MATERIAL DROP DOWN " --------------------------------------------------"""
#
# # Dropdown menu options
# options = [
#     "ABS - EU28+EFTA",
#     "ABS - non EU28+EFTA",
#     "Aluminum - Cu solute - EU28",
#     "Aluminum - Si solute - EU28",
#     "Aluminum - Mn solute - EU28+3",
#     "Aluminum - Mg solute - EU28+ETFA",
#     "Brass - EU28 +ETFA",
#     "Brass - Lead free - GLOBAL",
#     "Copper - EU28 +ETFA",
#     "Carbon Steel - Hot rolled coil - EU28+ETFA",
#     "Carbon Steel - Hot rolled coil - non EU28+ETFA",
#     "Carbon Steel - Cold rolled - ROW",
#     "Carbon Steel - Cold rolled coil - EU28+ETFA",
#     "NBR - non EU28+ETFA",
#     "NBR - non EU28+ETFA",
#     "PPS - EU28 +ETFA",
#     "PPS - non EU28 +ETFA",
#     "PP - EU28 +ETFA",
#     "PP - non EU28 +ETFA",
#     "PTFE - EU28 +ETFA",
#     "SS - Cold rolled - EU28+ETFA",
#     "SS - Hot rolled - ROW",
#     "SS - quatro plate - ROW",
#     "Zinc - GLOBAL",
#
# ]
#
# # datatype of menu text
# clicked = StringVar()
#
# # initial menu text
# clicked.set("Pick a material from the list")
# # Create Dropdown menu
#
# drop = OptionMenu(root, clicked, *options)
# drop.pack()
# # Create button, it will change label text
# button = Button(root, text="click Me", command=show).pack()
# # Create Label
# label = Label(root, text=" ")
# label.pack()
#
# # Execute tkinter
# root.mainloop()
# #
# """---------------------------------------------------------------------------"""



# """---------------------------------------------------------------------------"""


import tkinter as tk
from tkinter import filedialog

def UploadAction(event=None):
    filename = filedialog.askopenfilename()
    print('Selected:', filename)
    df_lca = pd.read_excel( filename, sheet_name="A1 - Raw material", index_col=0)

    return filename

root = tk.Tk()
button = tk.Button(root, text='Open', command=UploadAction)
button.pack()

root.mainloop()
# """---------------------------------------------------------------------------"""


#
# """ -----------Drop down then input as a separate box  ------------------ """
# import tkinter as tk
# from tkinter import simpledialog
#
#
# win = tk.Tk()
# win.geometry("100x50")
#
# def take_user_input_for_something():
#     user_input = simpledialog.askstring("Pop up for user input!", "What do you want to ask the user to input here?")
#     if user_input != "":
#         print(user_input)
#
# menubar = tk.Menu(win)
# dropDown = tk.Menu(menubar, tearoff = 0)
# dropDown.add_command(label = "Do something", command = take_user_input_for_something)
#
# # this entry field is not really needed her.
# # however I noticed you did not define this widget correctly
# # so I put this in to give you an example.
# my_entry = tk.Entry(win)
# my_entry.pack()
#
# menubar.add_cascade(label = "Drop Down", menu = dropDown)
#
#
#
# win.config(menu = menubar)
#
# win.mainloop()
# """---------------------------------------------------------------------------"""





# """---------------------------------------------------------------------------"""

# """---------------------------------------------------------------------------"""

# """---------------------------------------------------------------------------"""

# """---------------------------------------------------------------------------"""

# """---------------------------------------------------------------------------"""


from tkinter import *


class Table:

    def __init__(self, root):

        # code for creating table
        for i in range(total_rows):
            for j in range(total_columns):
                self.e = Entry(root, width=20, fg='black',
                               font=('Arial', 12, 'normal'))

                self.e.grid(row=i, column=j)
                self.e.insert(END, lst[i][j])


# take the data
lst = [ ("", 'Material type & source', 'Weight [gr]', "Average Recycling rate"),
        (1, 'material selected from dropdown menu ', 1, 19),
        (2, 'material selected from dropdown menu ', 2, 35),
        (3, 'material selected from dropdown menu ', 3, 35),
        (4, 'material selected from dropdown menu ', 4, 35),
        (5, 'material selected from dropdown menu ', 5, 35),
        ]


# find total number of rows and
# columns in list
total_rows = len(lst)
total_columns = len(lst[0])

# create root window
root = Tk()
t = Table(root)
root.mainloop()
# """---------------------------------------------------------------------------"""


