import paramiko
ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect('127.0.0.1', username='jesse', password='lol')


stdin, stdout, stderr = ssh.exec_command("uptime")

stdout.readlines()

def do_run(self, command):
        """run 
        Execute this command on all hosts in the list"""
        if command:
            for host, conn in zip(self.hosts, self.connections):
                stdin, stdout, stderr = conn.exec_command(command)
                stdin.close()
                for line in stdout.read().splitlines():
                    print 'host: %s: %s' % (host[0], line)
        else:
            print "usage: run "


def do_close(self, args):
    for conn in self.connections:
        conn.close()


if __name__ == '__main__':
    RunCommand().cmdloop()


# ftp = ssh.open_sftp() ftp.get('remotefile.py', 'localfile.py') ftp.close()
