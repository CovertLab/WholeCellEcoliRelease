
from django.conf.urls import patterns, url

from public import views

urlpatterns = patterns('',
	url(r'^$',views.index,name='index'),
	url(r'^list/(?P<model_type>\w+)/*$', views.list),
	url(r'^detail/(?P<model_type>\w+)/(?P<frame_id>[a-zA-Z0-9_\-]+)/*$', views.detail),
)

