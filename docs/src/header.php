<?php

require "password.php";

$passwd = explode("\n", file_get_contents('.htpasswd'));
$valid_passwords = array ();
foreach($passwd as $line)
{
    $line = preg_replace('/\s+/','',$line); // remove spaces
    if ($line) {
        list($user, $pass) = split(":", $line, 2);
        $valid_passwords[$user] = $pass;       
    }
}
$valid_users = array_keys($valid_passwords);

$user = $_SERVER['PHP_AUTH_USER'];
$pass = $_SERVER['PHP_AUTH_PW'];

$validated = (in_array($user, $valid_users)) && password_verify ($pass, $valid_passwords[$user]);

if (!$validated) {
  header('WWW-Authenticate: Basic realm="genome.crg.es"');
  header('HTTP/1.0 401 Unauthorized');
  die ("Not authorized");
}

?>
