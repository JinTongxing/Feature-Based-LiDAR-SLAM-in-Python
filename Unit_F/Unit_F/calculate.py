import string
import matplotlib.pyplot as plt

r = file('robot4_reference.txt')
r_data = r.readlines()
r_data_x = []
r_data_y = []
for rd in r_data:
    sp = rd.split()
    r_data_x.append(string.atof(sp[2]))
    r_data_y.append(string.atof(sp[3]))
r.close()

p = file('pure.txt')
p_data = p.readlines()
p_data_x = []
p_data_y = []
for pd in p_data:
    sp = pd.split()
    p_data_x.append(string.atof(sp[1]))
    p_data_y.append(string.atof(sp[2]))
p.close()

m = file('modify.txt')
m_data = m.readlines()
m_data_x = []
m_data_y = []
for md in m_data:
    sp = md.split()
    m_data_x.append(string.atof(sp[1]))
    m_data_y.append(string.atof(sp[2]))
m.close()

x = xrange(278)
y = []
z = []
# for i in xrange(278):
#     y.append(p_data_x[i] - r_data_x[i])
#     z.append(m_data_x[i] - np.array(r_data_x)
# y = abs(np.array(p_data_x) - np.array(r_data_x))
# z = abs(np.array(m_data_x) - np.array(r_data_x))
# p = abs(np.array(p_data_y) - np.array(r_data_y))
# q = abs(np.array(m_data_y) - np.array(r_data_y))
for i in xrange(278):
    y.append((p_data_x[i]-r_data_x[i])**2 + (p_data_y[i]-r_data_y[i])**2)
    z.append((m_data_x[i]-r_data_x[i])**2 + (m_data_y[i]-r_data_y[i])**2)
plt.figure()
plt.plot(x, y, 'r', linewidth=1)
plt.plot(x, z, 'b', linewidth=1)
plt.show()





