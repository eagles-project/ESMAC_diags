'''
colorblind color styles
'''


CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']


opacity = 0.5
colors = {
    'blue':   [55,  126, 184],  #377eb8 
    'orange': [255, 127, 0],    #ff7f00
    'green':  [77,  175, 74],   #4daf4a
    'pink':   [247, 129, 191],  #f781bf
    'brown':  [166, 86,  40],   #a65628
    'purple': [152, 78,  163],  #984ea3
    'gray':   [153, 153, 153],  #999999
    'red':    [228, 26,  28],   #e41a1c
    'yellow': [222, 222, 0]     #dede00
}  

#%%
# import matplotlib.pyplot as plt
# import numpy as np
# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1)
# x = np.array([0,1])
# for i in range(len(CB_color_cycle)):
#     slope = 1-i/9
#     ax.plot(x,x*slope,color=CB_color_cycle[i],label=i)
# ax.legend(loc='upper left')
