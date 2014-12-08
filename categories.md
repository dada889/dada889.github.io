layout: page
title: Categories
permalink: /Categories/
---
> 给我一双筷子，我能吃掉整个地球！ --------当代著名美食家**大衣申**
{% include JB/setup %}
<ul class="tag_box inline">
    {% assign categories_list = site.categories %}
    {% include JB/categories_list %}
</ul>
{% for category in site.categories %}
    <h2 id="{{ category[0] }}-ref">{{ category[0] | join: "/" }}</h2>
    <ul>
	{% assign pages_list = category[1] %}
	{% include JB/pages_list %}
    </ul>
{% endfor %}
