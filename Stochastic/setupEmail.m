function setupEmail()

    myaddress = 'mymatlabmail@gmail.com';
    mypassword = 'tp2fca8560_6988_422d_8745_3ec499db6668';

    setpref('Internet','E_mail',myaddress);
    setpref('Internet','SMTP_Server','smtp.gmail.com');
    setpref('Internet','SMTP_Username',myaddress);
    setpref('Internet','SMTP_Password',mypassword);

    props = java.lang.System.getProperties;
    props.setProperty('mail.smtp.auth','true');
    props.setProperty('mail.smtp.socketFactory.class', ...
                      'javax.net.ssl.SSLSocketFactory');
    props.setProperty('mail.smtp.socketFactory.port','465');

end

%sendmail(myaddress, 'Gmail Test', 'This is a test message.', 'emailTest.m');