# LPN-Estimator

## Introduction

This is an implementation of various algorithms designed to analyze and calculate the bit security of the Learning Parity with Noise (LPN) problem, as well as its dual problem Syndrome Decoding. The application supports analysis under exact or regular noise mode, and considers attacks including:

- Plain Gauss Elimination    
- Statistical Decoding and its [2.0](https://eprint.iacr.org/2022/1000.pdf)
- Information Set Decoding (including standard and BJMM variants)
- The algebraic attack in [this paper](https://eprint.iacr.org/2023/176)

It allows users to compute the complexity of various decoding strategies, providing insights about how secure their parameters are against several attacks.

Here is the parameters need for the LPN sample you want to analyze:

- $k$: the length of the secret
- $N$: the number of samples
- $t$: the hamming weight of the noise

And here are the parameters needed for the dual LPN, i.e., Syndrome Decoding sample you want to analyze:

- $n$: the number of correlations used in the dual LPN problem, another scale of the parity check matrix

You can also put your sample on larger fields for more algebraic structures and specify the field with:

- $q$: the size of the field if you want to use finite field $\mathbb{F}_q$
- $\lambda$: the logarithm size of the ring if you want to use ring $\mathbb Z_{2^\lambda}$

## Requirements

The application is written in Python 3.12.3 and requires the following packages:

```txt
asgiref==3.8.1
Django==5.1.4
gunicorn==23.0.0
numpy==2.2.1
packaging==24.2
sqlparse==0.5.3
```

which can be installed in a virtual environment using the following command:

```bash
pip install -r requirements.txt
```

## Deployment

If you want to simply run the django application without backend server, you can run the following command:

```bash
python manage.py runserver
```

If you want to deploy the application using Gunicorn and Nginx, you can follow the instructions below. Note that the static files have already be collected into the 'staticfiles' directory included in the pack, so you don't necessarily run `collectstatic` unless making some changes to the CSS or JavaScript files.

First navigate to the root directory and run the following command:

```bash
gunicorn --bind 0.0.0.0:8000 --timeout 120 lpnestimator.wsgi:application
```

or if you need different port, you can change the port number in the command above.

To configure Nginx to proxy requests to the application, add the following configuration to the server block in the Nginx configuration file (creating a file at `/etc/nginx/sites-available/lpnestimator`):

```nginx
server {
    listen 80;
    server_name your_domain_or_IP;

    location = /favicon.ico { access_log off; log_not_found off; }
    location /static/ {
        root /path/to/your/project;
    }

    location / {
        proxy_pass http://127.0.0.1:8000;  # Forward requests to Gunicorn
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }
}

```

and link the configuration file to the sites-enabled directory:

```bash
sudo ln -s /etc/nginx/sites-available/lpnestimator /etc/nginx/sites-enabled/
```

Test the Nginx configuration and
 restart Nginx to apply the changes:

```bash
sudo nginx -t
sudo systemctl restart nginx
```

## More to come
The transformation to C/C++ code is still in progress.