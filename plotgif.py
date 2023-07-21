import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

# 读取CSV文件
u = np.genfromtxt('u.csv', delimiter=',')

# 设置参数
fps = 30  # 帧速率
num_frames = u.shape[1]  # 总帧数

# 生成x轴坐标
x = np.linspace(0, 2, u.shape[0])

# 循环生成每一帧图像
frames = []
for i in range(int(num_frames/3)):
    # 绘制当前帧的图像
    plt.plot(x, u[:, 3*i])
    plt.xlabel('x')
    plt.ylabel('u')

    # 将当前帧的图像保存为PIL图像对象
    fig = plt.gcf()
    fig.canvas.draw()
    frame = Image.frombytes('RGB', fig.canvas.get_width_height(), fig.canvas.tostring_rgb())
    frames.append(frame)

    # 清除当前帧的图像
    plt.clf()

# 保存图像序列为GIF动画文件
frames[0].save('animation.gif', save_all=True, append_images=frames[1:], optimize=False, duration=int(1000 / fps), loop=0)
