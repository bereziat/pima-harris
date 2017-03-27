# blabla

harris: harris.c
	gcc $< -o $@ `imlib2-config --cflags --libs`  -lImlib2 -lm -L/opt/X11/lib -I/opt/X11/include
test1:
	./harris imgs/imagetest.png  0.001 
	open imageFinal.jpg

clean:
	rm -f imageFinal.jpg harris


