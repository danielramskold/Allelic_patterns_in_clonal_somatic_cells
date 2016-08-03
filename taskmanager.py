import threading, time, os, sys
from concurrent import futures

"""
Public classes and functions:
Tasklist(maxprocesses, verbose)
Tasklist.add(function, args, group, sample, waitfor, num_p, files_exist, maxingroup, kwargs)
Tasklist.waitforall()
Tasklist.anyerror()
Tasklist.groupresults(group)
Task.get()
"""

class Tasklist:
	def __init__(self, maxprocesses, verbose=False, allowoverrun=True, singleprocess=False):
		"""
		maxprocesses = how many processes are allowed at a time
		verbose = set to True if it should report progress to stdout
		allowoverrun = True if a process requesting N processes can run when M processes are free, 0<M<N
		singleprocess = True to disable multiprocessing and run functions in series
		"""		
		self.max_p = maxprocesses
		self.used_p = 0	# locked with lock_p
		self.tasks = []
		self.lock_p = threading.RLock()
		self.event_p = threading.Event()
		self.event_p.set()
		self.verbose = verbose
		self.printlock = threading.RLock()
		self.executor = futures.ProcessPoolExecutor(maxprocesses)
		self.allowoverrun = allowoverrun
		self.normal_msghandle = sys.stdout
		self.exception_msghandle = sys.stderr
		self.singleprocess = singleprocess
		self.abandonremaining = False
	
	def acquireprocesses(self, extra_p):
		bedtime = True
		while bedtime:
			with self.lock_p:
				if self.allowoverrun:
					bedtime = self.used_p != 0 and self.used_p >= self.max_p
				else:
					bedtime = self.used_p != 0 and self.used_p + extra_p > self.max_p
				if bedtime: self.event_p.clear()
				else: self.used_p += extra_p
			if self.abandonremaining:
				raise AbandonRemainingError
			if bedtime:
				self.event_p.wait()	# sleep until next releaseprocesses call
	
	def releaseprocesses(self, extra_p):
		with self.lock_p:
			self.used_p -= extra_p
			self.event_p.set()
	
	def waitforall(self):
		"""
		wait for all tasks to finish
		"""
		for task in self.tasks:
			task.thread.join()
			
	def anyerror(self):
		"""
		return True if any task has failed
		"""
		return any(task.badend for task in self.tasks)
	
	def groupresults(self, group):
		"""
		returns a generator of (sample, return value) tuples for tasks in given group
		"""
		return ((t.sample, t.get()) for t in self.tasks if t.group == group)
	
	def add(self, function, args=(), kwargs={}, group='unnamed', sample='all', waitfor=[], num_p=1, files_exist=[], maxingroup=None, retry_times=0, waitfor_moresamples=[], timeout=None):
		"""
		function = called function
		args = tuple of function arguments
		kwargs = dict of named function arguments
		any Task instance in args or kwargs will be substituted by Task.get() just before function call
		group = identifier (e.g. a name string) 
		sample = another identifier, set to 'all' to wait for all samples
		waitfor = list of groups that need to be run before, if they match the sample
		num_p = number of subprocesses spawned
		files_exist = list of files, cancel if any doesn't exists
		maxingroup = how many tasks can be in the same group before they wait for each other
		retry_times = how many times it will restart after an error (default: 0)
		timeout = how many seconds the function can run before it is force-crashed, ignored if tasklist.singleprocess == True
		returns a Task instance
		"""
		task = Task(function, args, self, group, sample, waitfor, num_p, files_exist, maxingroup, kwargs, retry_times, waitfor_moresamples, timeout)
		self.tasks.append(task)
		return task

	def sayevent(self, event, task):
		if self.verbose:
			with self.printlock:
				print >>self.normal_msghandle, ', '.join([event, str(task.group), str(task.sample), time.asctime()])

class Task:
	def __init__(self, function, args, tasklist, group, sample, waitfor, num_p, files_exist, maxingroup, kwargs, retry_times, waitfor_moresamples, timeout):
		waitfortasks = [task for task in tasklist.tasks if (sample == task.sample or sample=='all' or task.sample=='all' or task.sample in waitfor_moresamples) and task.group in waitfor]
		if maxingroup is not None:
			ingroup = [task for task in tasklist.tasks if task.group == group]
			if len(ingroup) >= maxingroup:
				waitfortasks.append(ingroup[-maxingroup])
		self.group = group
		self.sample = sample
		self.badend = False
		self.checkfiles = files_exist
		self.retries_left = retry_times + 1
		self.thread = threading.Thread(target=self.run, args=(function, args, kwargs, tasklist, num_p, waitfortasks, timeout))
		#tasklist.sayevent('Preparing', self)
		self.thread.start()
	
	def run(self, function, args, kwargs, tasklist, num_p, waitfortasks, timeout):
		for task in waitfortasks:
			task.thread.join()
			if task.badend:
				self.badend = True	# cascade
				tasklist.sayevent('Aborting', self)
				return
		try:
			# wait for Query instances' functions to finish and substitute with their return value
			args = tuple(arg.get() if isinstance(arg, Task) else arg for arg in args)
			kwargs = dict((key, arg.get()) if isinstance(arg, Task) else (key,arg) for key,arg in kwargs.items())
		except ValueError:
			self.badend = True
			tasklist.sayevent('Aborting', self)
			return
		for filename in self.checkfiles:
			if os.path.exists(filename): continue
			self.badend = True
			tasklist.sayevent('Aborting (missing ' + filename + ')', self)
			return
		completed = False
		while not completed and self.retries_left > 0:
			self.retries_left -= 1
			error_report = False, ''
			try: tasklist.acquireprocesses(num_p)
			except AbandonRemainingError:
				self.badend = True
				return
			try:
				tasklist.sayevent('Starting', self)
				time.sleep(0.01)
				if tasklist.singleprocess:
					error, ret = wrapperfunction(function, args, kwargs)
					if error:
						error_report = True, ret.rstrip()
					else:
						self.ret = ret
						completed = True
						tasklist.sayevent('Completed', self)
				elif function.__module__ != '__main__': # only wrap if the module name is known, otherwise there'll be an error
					pr = tasklist.executor.submit(wrapperfunction, function, args, kwargs)
					error, ret = pr.result(timeout=timeout)
					if error:
						error_report = True, ret.rstrip()
					else:
						self.ret = ret
						completed = True
						tasklist.sayevent('Completed', self)
				else:
					pr = tasklist.executor.submit(function, *args, **kwargs)
					try:
						self.ret = pr.result(timeout=timeout)
					except Exception as  e:
						error_report = True, e.__class__.__name__ + ': ' + str(e)
					else:
						completed = True
						tasklist.sayevent('Completed', self)
			except KeyboardInterruptError:
				tasklist.abandonremaining = True
				self.badend = True
				sys.stderr.flush()
				sys.stderr = open('/dev/null','w') # to skip the repeated tracebacks
			except KeyboardInterrupt:
				tasklist.abandonremaining = True
				self.badend = True
			finally:
				tasklist.releaseprocesses(num_p)
				if not completed:
					with tasklist.printlock:
						if error_report[0]:
							print >>tasklist.exception_msghandle, error_report[1]
						if self.retries_left:
							tasklist.sayevent('Error (tries left: %d)'%self.retries_left, self)
						else:
							tasklist.sayevent('Error', self)
							self.badend = True
				
	def get(self):
		"""
		waits for the task to finish and gives the function's return value
		"""		
		self.thread.join()
		if self.badend: raise ValueError
		return self.ret

def wrapperfunction(function, args, kwargs):
	try:
		r = function(*args, **kwargs)
		return 0, r
	except KeyboardInterrupt:
		raise KeyboardInterruptError
	except Exception:
		import traceback
		return 1, traceback.format_exc()

class KeyboardInterruptError(Exception): pass
class AbandonRemainingError(Exception): pass
