# Using `strkit visualize` remotely

The `strkit visualize` command starts a local web server in order to make genotyping results data browseable.
However, if you are working in a high-performance computing (HPC) environment such as a Digital Research Alliance of 
Canada (DRAC; f.k.a. Compute Canada) cluster, you will need to do some additional work to access this web interface 
without having to copy data to your local machine.

> [!NOTE]
> This guide assumes that the `visualize` web server will run on port 5011 both locally and on the server.

> [!IMPORTANT]
> While the server is running, others on the machine can access it, as it binds a port on the entire server.

As an example, let's assume we have some genotyping results data on the 
[Rorqual](https://docs.alliancecan.ca/wiki/Rorqual/en) cluster. We're going to forward the 5011 port to our local 
machine from whichever login node we sign in to, using the `ssh -L` option (with the syntax 
`<local port>:<remote host>:<remote port>`):

```bash
ssh -L 5011:localhost:5011 <user>@rorqual.alliancecan.ca 
cd path/to/our/results  # navigate to the directory with our genotyping results JSON file and other required files
source env/bin/activate  # activate our STRkit Python virtual environment
strkit visualize TODO
```

Then, on your local machine, go to http://localhost:5011. While the SSH connection is active and `strkit visualize` is 
running, you should be able to interact with the interface.
