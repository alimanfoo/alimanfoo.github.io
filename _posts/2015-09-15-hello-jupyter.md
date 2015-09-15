---
layout: post
title: Hello Jupyter
---


This is a test post written as a Jupyter notebook.


{% highlight python %}
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
{% endhighlight %}

Let's make a plot...


{% highlight python %}
x = np.random.random(1000)
plt.hist(x);
{% endhighlight %}


![png](/assets/2015-09-15-hello-jupyter_files/2015-09-15-hello-jupyter_3_0.png)


Let's make another plot...


{% highlight python %}
y = np.random.normal(size=1000)
plt.hist(y);
{% endhighlight %}


![png](/assets/2015-09-15-hello-jupyter_files/2015-09-15-hello-jupyter_5_0.png)


All done.
