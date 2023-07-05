FROM ubuntu

RUN apt-get clean && apt-get update && apt-get install -y -q cpanminus build-essential
RUN cpanm CGI
RUN cpanm Proc::Background
RUN cpanm File::Util
RUN cpanm File::Copy
RUN cpanm LWP::Simple

RUN apt update
RUN apt install -y apache2
RUN apt install -y apache2-utils
RUN apt clean

COPY ./000-default.conf /etc/apache2/sites-available/000-default.conf
COPY ./ports.conf /etc/apache2/ports.conf
COPY ./alias.conf /etc/apache2/mods-available/alias.conf

RUN a2enmod cgi.load
RUN a2enmod cgid.load

RUN mkdir -p /var/www/llama-magic
COPY . /var/www/llama-magic
RUN chown -hR www-data:www-data /var/www/llama-magic

EXPOSE 3000
CMD ["apache2ctl", "-D", "FOREGROUND"]


