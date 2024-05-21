from django.shortcuts import render
from django.http import HttpResponse
import home.estimator

import multiprocessing
import time
from multiprocessing import Process

cache = {}

def _calculator_with_timeout(request):
    

    keys = request.POST.copy()
    keys['csrfmiddlewaretoken']='none'

    if str(keys) in cache:
        return render(request, 'index.html', {'result': 'Load from Cache: \n' + cache[str(keys)]})

    manager = multiprocessing.Manager()
    res_dict = manager.dict()

    p = Process(target=_calculator, args=(request,res_dict))
    p.start()
    p.join(60)
    if p.is_alive():
        p.terminate()
        p.join()
        return render(request, 'index.html', {'result': "Timeout."})
    #return the output of process p
    

    cache[str(keys)] = res_dict['result']

    return  render(request, 'index.html', {'result': res_dict['result']})


def calculator(request):
    try:
        return _calculator_with_timeout(request)
    except Exception as e:
        print(e)
        return render(request, 'index.html', {'result': "An error occurred."})


def _calculator(request,res_dict):
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
    if N != None and N != "" and ((dual == 'on' and n != None and n != "") or (dual == None and k != None and k != "")) and t != None and t != "" and ((field == 'rq' and q != None and q != "") or (field == 'f2l' and l != None and l != "") or field == 'f2'):
        
        # convert to int
        N = int(N)
        if(dual == 'on'):
            n = int(n)
        else:
            k = int(k)
        t = int(t)
        if(field == 'rq'):
            q = int(q)
        if(field == 'f2l'):
            l = int(l)
            
        # compute the result
        if dual == 'on':
            if field == 'f2':
                if noise == 'exact':
                    result = "bit security of dual exact LPN (n=" + str(n) + ", N=" + str(N) + ", t=" + str(t) + "): " + str(home.estimator.analysisfordual2(n, N, t)) + " bits"
                elif noise == 'regular':
                    result = "bit security of dual regular LPN (n=" + str(n) + ", N=" + str(N) + ", t=" + str(t) + "): "+ str(home.estimator.analysisfordual2regular(n, N, t)) + " bits"
            elif field == 'rq':
                if noise == 'exact':
                    result = "bit security of dual exact LPN (n=" + str(n) + ", N=" + str(N) + ", t=" + str(t) + ", q=" + str(q) + "): " + str(home.estimator.analysisfordualq(n, N, t, q)) + " bits"
                elif noise == 'regular':
                    result = "bit security of regular LPN (n=" + str(n) + ", N=" + str(N) + ", t=" + str(t) + ", q=" + str(q) + "): " + str(home.estimator.analysisfordualqregular(n, N, t, q)) + " bits"
            elif field == 'f2l':
                if noise == 'exact':
                    result = "bit security of dual exact LPN (n=" + str(n) + ", N=" + str(N) + ", t=" + str(t) + ", lambda=" + str(
                    l) + "): " + str(home.estimator.analysisfordual2lambda(n, N, t, l)) + " bits"
                elif noise == 'regular':
                    result = "bit security of dual regular LPN (n=" + str(n) + ", N=" + str(N) + ", t=" + str(t) + ", lambda=" + str(
                    l) + "): " + str(home.estimator.analysisfordual2lambdaregular(n, N, t, l)) + " bits"
        else:
            if field == 'f2':
                if noise == 'exact':
                    result = "bit security of exact LPN (N=" + str(N) + ", k=" + str(k) + ", t=" + str(t) + "): " + str(home.estimator.analysisfor2(N, k, t)) + " bits"
                elif noise == 'regular':
                    result = "bit security of regular LPN (N=" + str(N) + ", k=" + str(k) + ", t=" + str(t) + "): " + str(home.estimator.analysisfor2regular(N, k, t)) + " bits"
            elif field == 'rq':
                if noise == 'exact':
                    result = "bit security of exact LPN (N=" + str(N) + ", k=" + str(k) + ", t=" + str(t) + ", q=" + str(q) + "): " + str(home.estimator.analysisforq(N, k, t, q)) + " bits"
                elif noise == 'regular':
                    result = "bit security of regular LPN (N=" + str(N) + ", k=" + str(k) + ", t=" + str(t) + ", q=" + str(q) + "): " + str(home.estimator.analysisforqregular(N, k, t, q)) + " bits"
            elif field == 'f2l':
                if noise == 'exact':
                    result = "bit security of exact LPN (N=" + str(N) + ", k=" + str(k) + ", t=" + str(t) + ", lambda=" + str(l) + "): " + str(home.estimator.analysisfor2lambda(N, k, t, l)) + " bits"
                elif noise == 'regular':
                    result = "bit security of regular LPN (N=" + str(N) + ", k=" + str(k) + ", t=" + str(t) + ", lambda=" + str(l) + "): " + str(home.estimator.analysisfor2lambdaregular(N, k, t, l)) + " bits"
        res_dict['result'] = result
        return render(request, 'index.html', {'result': result})
    else:
        res_dict['result'] = "Please fill out all fields."
        return render(request, 'index.html', {'result': "Please fill out all fields."})

def contact(request):
    return render(request, 'contact.html')