#!/bin/bash

set -x

base=$1
title=$2
ipynb=${base}.ipynb
md=${base}.md

echo Notebook: $ipynb
echo Title: $title
echo Output: $md 

# use nbconvert with vanilla markdown
rm -v $md
jupyter nbconvert --to markdown $ipynb

# insert frontmatter
echo -e "---\nlayout: post\ntitle: $title\n---\n" > ${md}.new
cat $md >> ${md}.new
mv ${md}.new ${md}

# fix up code blocks
sed -i -e 's/^```python$/{% highlight python %}/' $md
sed -i -e 's/^```$/{% endhighlight %}/' $md

