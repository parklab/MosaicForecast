#!/usr/bin/env python
import multiprocessing as mp
import time

fn = 'temp.txt'

def worker(arg, q):
    '''stupidly simulates long running process'''
    start = time.clock()
    s = 'this is a test'
    txt = s
    for i in range(1,2000):
        txt += s 
    done = time.clock() - start
    with open(fn, 'rb') as f:
        size = len(f.read())
    res = 'Process' + str(arg), str(size), done
    q.put(res)
    return res

def listener(q):
    '''listens for messages on the q, writes to file. '''

    f = open(fn, 'wb') 
    while 1:
        m = q.get()
        if m == 'kill':
            f.write('killed')
            break
        f.write(str(m) + '\n')
        f.flush()
    f.close()

def main():
    #must use Manager queue here, or will not work
    manager = mp.Manager()
    q = manager.Queue()    
    pool = mp.Pool(mp.cpu_count() + 2)

    #put listener to work first
    watcher = pool.apply_async(listener, (q,))

    #fire off workers
    jobs = []
    for i in range(80):
        job = pool.apply_async(worker, (i, q))
        jobs.append(job)

    # collect results from the workers through the pool result queue
    for job in jobs: 
        job.get()

    #now we are done, kill the listener
    q.put('kill')
    pool.close()

if __name__ == "__main__":
   main()
