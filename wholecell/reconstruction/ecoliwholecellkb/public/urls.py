
from django.conf.urls import patterns, url

from public import views

urlpatterns = patterns('',
	url(r'^$',views.index,name='index')
)

