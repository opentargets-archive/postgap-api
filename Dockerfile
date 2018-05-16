FROM node:carbon

# Create app directory
WORKDIR /usr/src/app

RUN yarn global add pm2

# Copy in the sqlite database
COPY postgap.20180514.db ./postgap.db
#TODO we should declare a VOLUME here

# Install app dependencies
# A wildcard is used to ensure both package.json AND package-lock.json are copied
# where available (npm@5+)
COPY package*.json ./

RUN yarn install
# If you are building your code for production
# RUN npm install --only=production

# Bundle app source
COPY . .

# Use production build
RUN yarn build

EXPOSE 4000
# CMD [ "yarn", "start" ]
ENTRYPOINT pm2 start --no-daemon dist/index.js