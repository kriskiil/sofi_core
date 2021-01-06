import { NextFunction, Request, Response } from 'express'
import config, { logger } from '../config'
import { CommonApi } from '@oryd/kratos-client'
import { IncomingMessage } from 'http'
import { isString } from '../helpers'

const kratos = new CommonApi(config.kratos.admin)

export default (req: Request, res: Response, next: NextFunction) => {
  const request = req.query.request

  // The request is used to identify the login and registration request and
  // return data like the csrf_token and so on.
  if (!request || !isString(request)) {
    logger.info('No request found in URL, initializing verify flow.')
    res.redirect(
      `${config.kratos.browser}/self-service/browser/flows/verification/email`
    )
    return
  }

  kratos
    .getSelfServiceVerificationRequest(request)
    .then(({ body, response }: { response: IncomingMessage; body?: any }) => {
      if (response.statusCode == 404) {
        res.redirect(
          `${config.kratos.browser}/self-service/browser/flows/verification/email`
        )
        return
      } else if (response.statusCode != 200) {
        return Promise.reject(body)
      }

      res.render('verification', body)
    })
    .catch(next)
}
