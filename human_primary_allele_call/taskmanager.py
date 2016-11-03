from __future__ import division, with_statement, print_function, unicode_literals, absolute_import
import threading, time, os, sys
from concurrent import futures

lastmodified = "7 Dec 2015"

"""
basic usage:
tasklist = taskmanager.Tasklist()
tasklist.add(your_function)
tasklist.add(your_function)
tasklist.run_all()

this is the taskmanager_lessthreaded version, built to work for both python2 and python3
"""

class Tasklist:
	def __init__(self, maxprocesses, verbose=False, singleprocess=False):
		"""
		maxprocesses = how many processes are allowed at a time
		verbose = set to True if it should report progress to stdout
		singleprocess = True to disable multiprocessing and run functions in series
		"""
		self.max_p = maxprocesses
		self.used_p = 0
		self.tasks = []
		self.verbose = verbose
		self.executor = futures.ProcessPoolExecutor(maxprocesses)
		self.normal_msghandle = sys.stdout
		#self.exception_msghandle = sys.stderr # not used
		self.singleprocess = singleprocess
		self.lock_p = threading.RLock()
		self.printlock = threading.RLock()
	
	def _allowed_tasks(self):
		for task_t in self.tasks:
			if task_t.done or task_t.running:
				continue
			if any((task in task_t.waitfor and not task.done) for task in self.tasks):
				continue
			if task_t.maxingroup is not None and sum(int(task.running and task.group == task_t.group) for task in self.tasks) >= task_t.maxingroup:
				continue
			if self.used_p + task_t.num_p * 0.7 > self.max_p: # allow minor overrun in terms of processes used, to not lock multi-process jobs behind single-process ones, thus the <1 coefficient
				continue
			yield task_t
	
	def run_all(self):
		"""
		runs all tasks added by add()
		"""
		times_waited = 0
		while not all(task.done for task in self.tasks):
			started_task = False
			with self.lock_p:
				if self.used_p < self.max_p:
					for task_t in self._allowed_tasks():
						if any(task.badend for task in task_t.waitfor):
							task_t.badend = True
							task_t.done = True
							self.sayevent('Abandoned', task_t)
						else:
							task_t._start()
							self.used_p += task_t.num_p
						started_task = True
						break
			if not started_task:
				time.sleep(0.05)
				times_waited += 1
				if times_waited == 10:
					if len(list(self._allowed_tasks())) == 0 and not any(task.running for task in self.tasks):
						self.sayevent('Waited for nothing')
						break
					times_waited = 0
			else:
				times_waited = 0
	
	def add(self, function, args=(), kwargs={}, group='unnamed', sample='all', waitfor=[], maxingroup=None, num_p=1, timeout_hours=None, waitfor_moresamples=[]):
		# missing argumens compared to taskmanager.py: files_exist, retry_times, timeout, require_wait
		"""
		adds a task, does not run it until run_all() is called
		
		function = called function
		args = tuple of function arguments
		kwargs = dict of named function arguments
		group = identifier (e.g. a name string) 
		sample = another identifier, set to 'all' to wait for all samples
		waitfor = list of groups that need to be run before, if they match the sample
		num_p = number of subprocesses spawned
		maxingroup = how many tasks can be in the same group before they wait for each other
		timeout_hours = how many hours the function can run before it is force-crashed, ignored if tasklist.singleprocess == True
		returns a Task instance
		"""
		task = Task(function, args, self, group, sample, waitfor, maxingroup, num_p, kwargs, timeout_hours, waitfor_moresamples)
		self.tasks.append(task)
		return task
		
	def sayevent(self, event, task=None):
		if self.verbose:
			with self.printlock:
				if task is None:
					print(', '.join([event, time.asctime()]), file=self.normal_msghandle)
				else:
					print(', '.join([event, str(task.group), str(task.sample), time.asctime()]), file=self.normal_msghandle)
	
	def waitforall(self):
		self.run_all()
	
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

class Task:
	def __init__(self, function, args, tasklist, group, sample, waitfor, maxingroup, num_p, kwargs, timeout_hours, waitfor_moresamples):
		"""
		called by running Tasklist.add(), not meant for calling directly
		"""
		self.function = function
		self.group = group
		self.sample = sample
		self.badend = False
		self.running = False
		self.done = False
		self.args = args
		self.kwargs = kwargs
		self.maxingroup = maxingroup
		self.tasklist = tasklist
		self.num_p = num_p
		self.timeout = None if timeout_hours is None else 3600*timeout_hours
		self.waitfor = set(task for task in tasklist.tasks if (sample == task.sample or sample=='all' or task.sample=='all' or task.sample in waitfor_moresamples) and task.group in waitfor)
	
	def _start(self):
		self.running = True
		thread = threading.Thread(target=self._run_in_thread, args=(self.function, self.args, self.kwargs, self.tasklist))
		thread.start()
	
	def _run_in_thread(self, function, args, kwargs, tasklist):
		tasklist.sayevent('Starting', self)
		try:
			if tasklist.singleprocess:
				try:
					self.ret = function(*args, **kwargs)
				except:
					self.badend = True
					tasklist.sayevent('Error', self)
					raise
			else:
				pr = tasklist.executor.submit(function, *args, **kwargs)
				try:
					self.ret = pr.result(timeout=self.timeout)
				except Exception as e:
					self.badend = True
					self.tasklist.sayevent('Error: %s'%repr(e), self)
			if not self.badend:
				tasklist.sayevent('Completed', self)
		finally:
			with tasklist.lock_p:
				self.running = False
				self.done = True
				self.tasklist.used_p -= self.num_p
	
	def get(self):
		"""
		waits for the task to finish and gives the function's return value
		"""
		if not self.done:
			self.tasklist.run_all()
			if self.badend: raise ValueError
		return self.ret
