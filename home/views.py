from django.shortcuts import render
from django.http import HttpResponse
import home.estimator
 
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

    result = ""

    # convert to int
    if N == '' or N == None:
        N = 0
    else:
        N = int(N)
    if n == '' or n == None:
        n = 0
    else:
        n = int(n)
    if k == '' or k == None:
        k = 0
    else:
        k = int(k)
    if t == '' or t == None:
        t = 0
    else:
        t = int(t)
    if q == '' or q == None:
        q = 0
    else:
        q = int(q)
    if l == '' or l == None:
        l = 0
    else:
        l = int(l)
        
    # compute the result
    if dual == 'on':
        if field == 'f2':
            if noise == 'exact':
                result = "bit security of dual exact LPN (n=" + str(n) + ", N=" + str(N) + ", t=" + str(t) + "):" + str(home.estimator.analysisfordual2(n, N, t)) + " bits"
            elif noise == 'regular':
                result = "bit security of dual regular LPN (n=" + str(n) + ", N=" + str(N) + ", t=" + str(t) + "):"+ str(home.estimator.analysisfordual2regular(n, N, t)) + " bits"
        elif field == 'rq':
            if noise == 'exact':
                result = "bit security of dual exact LPN (n=" + str(n) + ", N=" + str(N) + ", t=" + str(t) + ", q=" + str(q) + "):" + str(home.estimator.analysisfordualq(n, N, t, q)) + " bits"
            elif noise == 'regular':
                result = "bit security of regular LPN (n=" + str(n) + ", N=" + str(N) + ", t=" + str(t) + ", q=" + str(q) + "):" + str(home.estimator.analysisfordualqregular(n, N, t, q)) + " bits"
        elif field == 'f2l':
            if noise == 'exact':
                result = "bit security of dual exact LPN (n=" + str(n) + ", N=" + str(N) + ", t=" + str(t) + ", lambda=" + str(
                l) + "):" + str(home.estimator.analysisfordual2lambda(n, N, t, l)) + " bits"
            elif noise == 'regular':
                result = "bit security of dual regular LPN (n=" + str(n) + ", N=" + str(N) + ", t=" + str(t) + ", lambda=" + str(
                l) + "):" + str(home.estimator.analysisfordual2lambdaregular(n, N, t, l)) + " bits"
    else:
        if field == 'f2':
            if noise == 'exact':
                result = "bit security of exact LPN (N=" + str(N) + ", k=" + str(k) + ", t=" + str(t) + "):" + str(home.estimator.analysisfor2(N, k, t)) + " bits"
            elif noise == 'regular':
                result = "bit security of regular LPN (N=" + str(N) + ", k=" + str(k) + ", t=" + str(t) + "):" + str(home.estimator.analysisfor2regular(N, k, t)) + " bits"
        elif field == 'rq':
            if noise == 'exact':
                result = "bit security of exact LPN (N=" + str(N) + ", k=" + str(k) + ", t=" + str(t) + ", q=" + str(q) + "):" + str(home.estimator.analysisforq(N, k, t, q)) + " bits"
            elif noise == 'regular':
                result = "bit security of regular LPN (N=" + str(N) + ", k=" + str(k) + ", t=" + str(t) + ", q=" + str(q) + "):" + str(home.estimator.analysisforqregular(N, k, t, q)) + " bits"
        elif field == 'f2l':
            if noise == 'exact':
                result = "bit security of exact LPN (N=" + str(N) + ", k=" + str(k) + ", t=" + str(t) + ", lambda=" + str(l) + "):" + str(home.estimator.analysisfor2lambda(N, k, t, l)) + " bits"
            elif noise == 'regular':
                result = "bit security of regular LPN (N=" + str(N) + ", k=" + str(k) + ", t=" + str(t) + ", lambda=" + str(l) + "):" + str(home.estimator.analysisfor2lambdaregular(N, k, t, l)) + " bits"
    return render(request, 'index.html', {'result': result})