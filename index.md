---
layout: default
---

> Talk is cheap. Show me the code. ----- Linus Torvalds


<div class="home">
<!--
  <h1 class="page-heading">My Articles</h1>
-->

  <ul class="post-list">
    {% for post in site.posts %}
      <li>
        <span class="post-meta">{{ post.date | date: "%b %-d, %Y" }}</span>

        <h2>
          <a class="post-link" href="{{ post.url | prepend: site.baseurl }}">{{ post.title }}</a>
        </h2>
	<hr style="color:#EEEEFF" />
	<p>{{post.content}}</p>
	
      </li>
    {% endfor %}
  </ul>

  <p class="rss-subscribe">subscribe <a href="{{ "/feed.xml" | prepend: site.baseurl }}">via RSS</a></p>

</div>
{% include footer.html %}
