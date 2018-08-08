import json
from confluent_kafka import Producer, Consumer, KafkaError

def delivery_report(err, msg):
	"""
	This is a utility method passed to the Kafka Producer to handle the delivery
	of messages sent using `send(topic, message)`. 
	"""
	if err is not None:
		print('message delivery failed: {}'.format(msg))
		print('failed message: {}'.format(err))

class Agent(object):

	"""
	Agent: a process with the ability to send and receive messages using Kafka.

	This is a base class that handles all of the details of connecting to both a 
	Kafka Producer (to send messages) and Consumer (to receive messages). 

	To subclass Agent:

	* First, a call to `super(Subclass, self).__init__(id, kafka)` must be made in the
      `__init__` method of the overriding class. `id` is any unique string and identifies
      this agent with other agents in the system. `kafka` is a configuration dictionary
      containing information about the kafka host, topics to subscribe to and additional 
      information. Its details are described in the docstring for `__init__`.

	* Then, you can override `initialize()` to send any messages before the consumer
	  loop begins.

	* The main method to override is `receive(topic, message)`, which will be called
	  each time a message is received on any of the topics originally subscribed to.

	* Finally, on shutdown there is a call to `finalize()` where you can put anything
	  that needs to be cleaned up or any final messages to send before the process is
	  terminated.

	Additionally, a `send(topic, message)` method is provided that does not need to be
	overridden but can be used directly, which sends a message to the specified topic
	using the configured Kafka Producer.
	"""

	def __init__(self, id, kafka):
		"""
		Initialize the Agent object with a unique id and kafka configuration.

		Args:
		    id (str): A unique identifier which distinguishes this agent from 
		        the rest of the agents in the system.

		    kafka (dict): A dictionary containing Kafka configuration information.
		        The relevant keys are:

		        * host (str): The address of the Kafka host.
		        * subscribe_topics (array[str]): A list of topics to subscribe to. 
		            If this is an empty array then no consumer will be initialized.

                Additional keys can be provided for use in the overriding agent class
		        to represent topics for sending messages to. 
		"""

		self.id = id
		self.kafka = kafka

		self.producer = Producer({
			'bootstrap.servers': self.kafka['host']})

		self.has_consumer = len(self.kafka['subscribe_topics']) > 0
		if self.has_consumer:
			self.consumer = Consumer({
				'bootstrap.servers': self.kafka['host'],
				'enable.auto.commit': True,
				'group.id': 'simulation-' + str(id),
				'default.topic.config': {
					'auto.offset.reset': 'latest'}})

		self.initialize()

		if self.has_consumer:
			self.consumer.subscribe(
				self.kafka['subscribe_topics'])

			self.poll()

	def initialize(self):
		"""
		Initialize the Agent in the system.

		This method is called after the Producer has been initialized but before the 
		consumer loop is started. It is a good place to send any initialization messages
		to other agents before falling into the polling cycle.
		"""

		pass

	def poll(self):
		"""
		Enter the main consumer polling loop.

		Once poll is called, the thread will be claimed and any interaction with the 
		system from this point on will be mediated through message passing. This is called
		at the end of the base class's `__init__(id, kafka)` method and does not need to be
		called manually by the subclass.
		"""

		self.running = True
		while self.running:
			raw = self.consumer.poll()
			if raw is None:
				continue
			if raw.error():
				if raw.error().code() == KafkaError._PARTITION_EOF:
					continue
				else:
					print(raw.error())
					self.running = False

			message = json.loads(raw.value().decode('utf-8'))

			if message['event'] == 'GLOBAL_SHUTDOWN':
				self.shutdown()
			else:
				self.receive(raw.topic(), message)

	def send(self, topic, message):
		"""
		Send a dictionary as a message on the given topic. 

		Args:
		    topic (str): The Kafka topic to send the message on.
		    message (dict): A dictionary containing the message to send. This dictionary
		        needs to be JSON serializable, so must contain only basic types like `str`,
		        `int`, `float`, `list`, `tuple`, `array`, and `dict`. Any functions or objects
                present will throw errors.
		"""

		print('<-- ' + topic + ': ' + str(message))

		self.producer.flush()
		self.producer.poll(0)
		self.producer.produce(
			topic,
			json.dumps(message).encode('utf-8'),
			callback=delivery_report)

	def receive(self, topic, message):
		"""
		Receive a message from the given topic.

		This is the main method that needs to be overridden by a subclass to provide the 
		behavior for the agent. This method is invoked each time the agent receives a message
		on one of the topics originally subscribed to, which was supplied in the 
		kafka['subscribe_topics'] array.

		By convention, each message contains an `event` key which can be switched on in this
		method to trigger a specific response.

		Args:
		    topic (str): The Kafka topic this message was received on.
		    message (dict): A deserialized dictionary containing at the least an `event` key
		        describing the event that originally triggered this message.
		"""

		pass

	def shutdown(self):
		"""
		Shutdown the Kafka Producer and Consumer.

		Breaks out of the polling receive loop and calls `finalize`, then flushes the
		producer, commits the consumer offsets and exits.
		"""

		self.running = False
		self.finalize()

		self.producer.flush()
		self.producer.poll(0)

		if self.has_consumer:
			self.consumer.commit()
			self.consumer.close()

	def finalize(self):
		"""
		Send any final messages and do any clean up before exiting.

		This method can be overridden by the subclass to send any final messages to other 
		agents in the system or release system resources before the agent is shut down.
		"""

		pass
