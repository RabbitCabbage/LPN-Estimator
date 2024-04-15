from django.shortcuts import render
from django.http import HttpResponse
 
def calculator(request):
    noise = request.POST.get('noise') # exact or regular
    dual = request.POST.get('dual') # on or None
    field = request.POST.get('field') # f2 rq or f2l
    N, n, k, t, q, l = 0, 0, 0, 0, 0, 0
    N = request.POST.get('N')
    if(dual == 'on'):
        n = request.POST.get('n')
    k = request.POST.get('k')
    t = request.POST.get('t')
    if(field == 'rq'):
        q = request.POST.get('q')
    if(field == 'f2l'):
        l = request.POST.get('l')
    return render(request, 'index.html', {'result': N})