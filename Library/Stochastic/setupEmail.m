function setupEmail()

    myaddress = 'yourAddressGoesHere';
    mypassword = 'yourPasswordGoesHere';
    
    % this may vary if you're not using gmail
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

% test sending email
%sendmail(myaddress, 'Gmail Test', 'This is a test message.', 'emailTest.m');