from tkinter import *
from PIL import Image,ImageTk

root = Tk()

path = '/Users/muyuyang/Desktop/Frith\ Lab/results/oryza\ chloroplast\ Sativa\ vs\ Coarctata/os-oc-2.png'.replace('\\','')

image = Image.open(path)

# image = image.resize((500,500))
newImage = ImageTk.PhotoImage(image)


canvas = Canvas(root,width=1000,height=1000)
canvas.pack()
canvas.create_image(0,0,image=newImage)
canvas.create_rectangle(0,0,500,500,outline = 'blue')

root.mainloop()