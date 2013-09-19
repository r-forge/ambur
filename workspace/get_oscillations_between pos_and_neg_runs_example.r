x = c(-23 , -2 , 5 , 8, 9, 12, -2, -1, 3, 5, 7)

length(rle(sign(x))$lengths)  # gets each oscillation from between positive and negative

length(x) #get the number of eras

(rle(sign(x))$values) #get the positive vs negative oscillation


sum((rle(sign(x))$values * rle(sign(x))$lengths)[(rle(sign(x))$values * rle(sign(x))$lengths)<0]) #get number of eras with erosion
